
"""Cell ranger scRNA-Seq analysis."""
import os

from shutil import move
from resolwe_runtime_utils import export
from resolwe.process import (
    Cmd,
    DataField,
    DirField,
    FileField,
    FileHtmlField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
    UrlField,
)


class CellRangerMkref(Process):
    """Reference preparation tool for 10x Genomics Cell Ranger.

    Build a Cell Ranger-compatible reference from genome FASTA and gene GTF files.
    https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

    """

    slug = 'cellranger-mkref'
    name = 'Cell Ranger Mkref'
    process_type = 'data:genomeindex:10x'
    version = '1.0.0'
    category = 'scRNA-Seq'
    sheduling_class = SchedulingClass.BATCH
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/scseq:1.0.0'
            }
        },
        'resources': {
            'memory': 32768,
            'cores': 10,
        },
    }
    data_name = '{{ genome.fasta.file|default("?") }}'

    class Input:
        """Input fields to process CellRangerMkref."""

        genome = DataField(
            data_type='genome:fasta:',

            label='Reference genome',
        )
        annotation = DataField(
            data_type='annotation:gtf:',
            label='Annotation',
        )

    class Output:
        """Output fields to process CellRangerMkref."""

        genome_index = DirField(label="Indexed genome")
        build = StringField(label="Build")
        species = StringField(label="Species")
        source = StringField(label="Gene ID source")

    def run(self, inputs, outputs):
        """Run the analysis."""
        genome_build = inputs.genome.build
        annotation_build = inputs.annotation.build
        if genome_build != annotation_build:
            self.error(
                'Builds of the genome {} and annotation {} do not match. Please provide genome '
                'and annotation with the same build.'.format(genome_build, annotation_build)
            )

        genome_species = inputs.genome.species
        annotation_species = inputs.annotation.species
        if genome_species != annotation_species:
            self.error(
                'Species of genome {} and annotation {} do not match. Please provide genome '
                'and annotation with the same species.'.format(genome_species, annotation_species)
            )

        cmd = Cmd['cellranger']['mkref']
        cmd = cmd['--genome={}'.format(genome_build)]
        cmd = cmd['--genes={}'.format(inputs.annotation.annot_sorted.path)]
        cmd = cmd['--fasta={}'.format(inputs.genome.fasta.path)]
        cmd = cmd['--nthreads={}'.format(self.requirements.resources.cores)]
        cmd = cmd['--memgb={}'.format(int(self.requirements.resources.memory * 0.9 / 1024))]
        cmd()

        os.rename(genome_build, 'cellranger_index')

        outputs.genome_index = 'cellranger_index'
        outputs.source = inputs.annotation.source
        outputs.species = genome_species
        outputs.build = genome_build


class CellRangerCount(Process):
    """Perform gene expression analysis.

    Generate single cell feature counts for a single library.
    https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

    """

    slug = 'cellranger-count'
    name = 'Cell Ranger Count'
    process_type = 'data:scexpression:10x'
    version = '1.0.0'
    category = 'scRNA-Seq'
    sheduling_class = SchedulingClass.BATCH
    entity = {
        'type': 'sample'
    }
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/scseq:1.0.0'
            }
        },
        'resources': {
            'memory': 32768,
            'cores': 10,
        },
    }
    data_name = '{{ reads|sample_name|default("?") }}'

    class Input:
        """Input fields to process ImportScRNA10x."""

        reads = DataField(
            data_type='screads:10x:',
            label='10x reads data object',
        )
        genome_index = DataField(
            data_type='genomeindex:10x:',
            label='10x genome index data object',
        )
        chemistry = StringField(
            label='Chemistry',
            required=False,
            default='auto',
            description=('Assay configuration. By default the assay configuration is detected '
                         'automatically, which is the recommended mode. You should only specify '
                         'chemistry if there is an error in automatic detection.'),
            choices=[
                ("auto", "auto"),
                ("Single Cell 3′", "threeprime"),
                ("Single Cell 5′", "fiveprime"),
                ("Single Cell 3′ v1", "SC3Pv1"),
                ("Single Cell 3′ v2", "SC3Pv2"),
                ("Single Cell 3′ v3", "SC3Pv3"),
                ("Single Cell 5′ paired-end", "C5P-PE"),
                ("Single Cell 5′ R2-only", "SC5P-R2"),
            ],
        )
        trim_r1 = IntegerField(
            label='Trim R1',
            required=False,
            description=('Hard-trim the input R1 sequence to this length. Note that the length '
                         'includes the Barcode and UMI sequences so do not set this below 26 for '
                         'Single Cell 3′ v2 or Single Cell 5′. This and "Trim R2" are useful for '
                         'determining the optimal read length for sequencing.'),
        )
        trim_r2 = IntegerField(
            label='Trim R2',
            required=False,
            description='Hard-trim the input R2 sequence to this length.',
        )
        expected_cells = IntegerField(
            label='Expected number of recovered cells',
            default=3000,
        )
        force_cells = IntegerField(
            label='Force cell number',
            required=False,
            description=('Force pipeline to use this number of cells, bypassing the cell '
                         'detection algorithm. Use this if the number of cells estimated by Cell '
                         'Ranger is not consistent with the barcode rank plot.'),
        )

    class Output:
        """Output fields to process ImportScRNA10x."""

        matrix_filtered = FileField(label='Matrix (filtered)')
        genes_filtered = FileField(label='Genes (filtered)')
        barcodes_filtered = FileField(label='Barcodes (filtered)')
        matrix_raw = FileField(label='Matrix (raw)')
        genes_raw = FileField(label='Genes (raw)')
        barcodes_raw = FileField(label='Barcodes (raw)')
        report = FileHtmlField(label="Report")
        build = StringField(label="Build")
        species = StringField(label="Species")
        source = StringField(label="Gene ID source")

    def run(self, inputs, outputs):
        """Run the analysis."""
        sample_name = inputs.reads.entity_name.strip('.fastq.gz')

        dir_fastqs = './fastqs'
        os.mkdir(dir_fastqs)

        # Format cellranger count fastq input so it follows the correct naming convention and
        # folder structure
        for i, fastqs in enumerate(zip(inputs.reads.barcodes, inputs.reads.reads)):
            os.symlink(
                fastqs[0].path,
                os.path.join(
                    dir_fastqs,
                    '{}_S1_L{}_R1_001.fastq.gz'.format(
                        sample_name,
                        str(i + 1).zfill(3),
                    )
                )
            )
            os.symlink(
                fastqs[1].path,
                os.path.join(
                    dir_fastqs,
                    '{}_S1_L{}_R2_001.fastq.gz'.format(
                        sample_name,
                        str(i + 1).zfill(3),
                    )
                )
            )

        cmd = Cmd['cellranger']['count']
        cmd = cmd['--id={}'.format(sample_name)]
        cmd = cmd['--fastqs={}'.format(dir_fastqs)]
        cmd = cmd['--transcriptome={}'.format(inputs.genome_index.genome_index.path)]
        cmd = cmd['--localcores={}'.format(self.requirements.resources.cores)]
        cmd = cmd['--localmem={}'.format(int(self.requirements.resources.memory * 0.9 / 1024))]
        cmd = cmd['--chemistry={}'.format(inputs.chemistry)]
        cmd = cmd['--expect-cells={}'.format(inputs.expected_cells)]
        if inputs.trim_r1:
            cmd = cmd['--r1-length={}'.format(inputs.trim_r1)]
        if inputs.trim_r2:
            cmd = cmd['--r2-length={}'.format(inputs.trim_r2)]
        if inputs.force_cells:
            cmd = cmd['--force-cells={}'.format(inputs.force_cells)]
        cmd()

        output_dir = '{}/outs'.format(sample_name)

        # Spawn upload-bam process
        bam_name = '{}.bam'.format(sample_name)
        bai_name = '{}.bam.bai'.format(sample_name)
        bam_path = os.path.join(output_dir, 'possorted_genome_bam.bam')
        bai_path = os.path.join(output_dir, 'possorted_genome_bam.bam.bai')
        print(export(move(bam_path, bam_name)))  # TODO: fix when self.export() is available
        print(export(move(bai_path, bai_name)))
        process_inputs = {
            'src': bam_name,
            'src2': bai_name,
            'reads': inputs.reads.id,
            'species': inputs.genome_index.species,
            'build': inputs.genome_index.build,
        }
        self.run_process('upload-bam-scseq-indexed', process_inputs)

        filtered_dir = os.path.join(output_dir, 'filtered_feature_bc_matrix')
        raw_dir = os.path.join(output_dir, 'raw_feature_bc_matrix')
        outputs.matrix_filtered = os.path.join(filtered_dir, 'matrix.mtx.gz')
        outputs.genes_filtered = os.path.join(filtered_dir, 'features.tsv.gz')
        outputs.barcodes_filtered = os.path.join(filtered_dir, 'barcodes.tsv.gz')
        outputs.matrix_raw = os.path.join(raw_dir, 'matrix.mtx.gz')
        outputs.genes_raw = os.path.join(raw_dir, 'features.tsv.gz')
        outputs.barcodes_raw = os.path.join(raw_dir, 'barcodes.tsv.gz')
        outputs.report = os.path.join(output_dir, 'web_summary.html')
        outputs.build = inputs.genome_index.build
        outputs.species = inputs.genome_index.species
        outputs.source = inputs.genome_index.source

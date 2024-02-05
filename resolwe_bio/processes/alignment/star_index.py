"""Generate STAR genome index."""

import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    DirField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    StringField,
)


class StarIndex(Process):
    """Generate STAR genome index.

    Generate genome indices files from the supplied reference genome
    sequence and GTF files. The current version of STAR is 2.7.10b.
    """

    slug = "alignment-star-index"
    name = "STAR genome index"
    process_type = "data:index:star"
    version = "4.0.0"
    category = "Genome index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.2.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 32768,
        },
    }
    data_name = "{{ ref_seq.fasta.file|basename|default('?') }}"

    class Input:
        """Input fields to process StarIndex."""

        ref_seq = DataField(
            "seq:nucleotide", label="Reference sequence (nucleotide FASTA)"
        )
        annotation = DataField(
            "annotation",
            label="Annotation file (GTF/GFF3)",
            required=False,
            description="Insert known annotations into genome indices at the indexing stage.",
        )
        source = StringField(
            label="Gene ID Database Source",
            disabled="annotation",
            required=False,
            allow_custom_choice=True,
            choices=[
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )

        class AnnotationOptions:
            """Annotation file options."""

            feature_exon = StringField(
                label="Feature type [--sjdbGTFfeatureExon]",
                default="exon",
                description="Feature type in GTF file to be used as exons for building "
                "transcripts.",
            )
            sjdb_overhang = IntegerField(
                label="Junction length [--sjdbOverhang]",
                default=100,
                description="This parameter specifies the length of the genomic sequence around "
                "the annotated junction to be used in constructing the splice junction database. "
                "Ideally, this length should be equal to the ReadLength-1, where ReadLength is "
                "the length of the reads. For instance, for Illumina 2x100b paired-end reads, the "
                "ideal value is 100-1=99. In case of reads of varying length, the ideal value is "
                "max(ReadLength)-1. In most cases, the default value of 100 will work as well as "
                "the ideal value.",
            )

        class AdvancedOptions:
            """Advanced options."""

            genome_sa_string_len = IntegerField(
                label="Small genome adjustment [--genomeSAindexNbases]",
                required=False,
                description="For small genomes, the parameter --genomeSAindexNbases needs to be "
                "scaled down, with a typical value of min(14, log2(GenomeLength)/2 - 1). For "
                "example, for 1 megaBase genome, this is equal to 9, for 100 kiloBase genome, "
                "this is equal to 7.",
            )
            genome_chr_bin_size = IntegerField(
                label="Bin size for genome storage [--genomeChrBinNbits]",
                required=False,
                description="If you are using a genome with a large (>5,000) number of references "
                "(chrosomes/scaffolds), you may need to reduce the --genomeChrBinNbits to reduce "
                "RAM consumption. The following scaling is recommended: --genomeChrBinNbits = "
                "min(18, log2(GenomeLength / NumberOfReferences)). For example, for 3 gigaBase "
                "genome with 100,000 chromosomes/scaffolds, this is equal to 15.",
            )

            genome_sa_sparsity = IntegerField(
                label="Suffix array sparsity [--genomeSAsparseD]",
                required=False,
                description="Suffix array sparsity, i.e. distance between indices: use bigger "
                "numbers to decrease needed RAM at the cost of mapping speed reduction (integer > "
                "0, default = 1).",
            )

        annotation_options = GroupField(
            AnnotationOptions, label="Annotation file options", hidden="!annotation"
        )
        advanced = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields to process StarIndex."""

        index = DirField(label="Indexed genome")
        fastagz = FileField(label="FASTA file (compressed)")
        fasta = FileField(label="FASTA file")
        fai = FileField(label="FASTA file index")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        if not inputs.source and not inputs.annotation:
            self.error(
                "Gene ID database Source information must be provided when annotation GTF is not "
                "selected."
            )

        fasta = Path(inputs.ref_seq.output.fasta.path).name
        shutil.copy(inputs.ref_seq.output.fasta.path, fasta)
        outputs.fasta = fasta

        fastagz = Path(inputs.ref_seq.output.fastagz.path).name
        shutil.copy(inputs.ref_seq.output.fastagz.path, fastagz)
        outputs.fastagz = fastagz

        fai = Path(inputs.ref_seq.output.fai.path).name
        shutil.copy(inputs.ref_seq.output.fai.path, fai)
        outputs.fai = fai

        self.progress(0.1)

        index_dir = Path("star_index")
        index_dir.mkdir()

        index_params = [
            "--runThreadN",
            self.requirements.resources.cores,
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(index_dir),
            "--genomeFastaFiles",
            str(Path(inputs.ref_seq.output.fasta.path).name),
        ]

        if inputs.annotation:
            index_params.extend(
                [
                    "--sjdbGTFfile",
                    inputs.annotation.output.annot.path,
                    "--sjdbOverhang",
                    inputs.annotation_options.sjdb_overhang,
                    "--sjdbGTFfeatureExon",
                    inputs.annotation_options.feature_exon,
                ]
            )

            if inputs.annotation.type.startswith("data:annotation:gff3:"):
                index_params.extend(["--sjdbGTFtagExonParentTranscript Parent"])

        if inputs.advanced.genome_sa_string_len:
            index_params.extend(
                ["--genomeSAindexNbases", inputs.advanced.genome_sa_string_len]
            )

        if inputs.advanced.genome_chr_bin_size:
            index_params.extend(
                ["--genomeChrBinNbits", inputs.advanced.genome_chr_bin_size]
            )

        if inputs.advanced.genome_sa_sparsity:
            index_params.extend(
                ["--genomeSAsparseD", inputs.advanced.genome_sa_sparsity]
            )

        return_code, _, _ = Cmd["STAR"][index_params] & TEE(retcode=None)
        if return_code:
            self.error("Genome index build failed.")

        self.progress(0.8)

        outputs.index = str(index_dir)

        if inputs.annotation:
            outputs.source = inputs.annotation.output.source
            outputs.species = inputs.annotation.output.species
            outputs.build = inputs.annotation.output.build
        else:
            outputs.source = inputs.source
            outputs.species = inputs.ref_seq.output.species
            outputs.build = inputs.ref_seq.output.build

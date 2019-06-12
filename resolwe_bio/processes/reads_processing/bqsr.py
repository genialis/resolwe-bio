"""Base quality score recalibration."""
import os
import shutil

from resolwe.process import Process, Cmd, SchedulingClass, DataField, \
    FileField, StringField, ListField


class BQSR(Process):
    """
    A two pass process of BaseRecalibrator and ApplyBQSR from GATK.

    See [GATK website](https://software.broadinstitute.org/gatk/documentation/tooldocs/current
    /org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) for more information.

    It is possible to modify read group using [AddOrReplaceGroups](
    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.11.0/picard_sam_AddOrReplaceReadGroups.php)
    through ``read_group`` input field.
    """

    slug = 'bqsr'
    name = 'BaseQualityScoreRecalibrator'
    process_type = 'data:alignment:bam:bqsr:'
    version = '1.0.0'
    category = 'BAM processing'
    shaduling_class = SchedulingClass.BATCH
    entity = {'type': 'sample'}
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/dnaseq:4.2.0'
            }
        },
    }
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields to perform Base quality score recalibration."""

        bam = DataField('alignment:bam', label='BAM file containing reads')
        reference = DataField('genome:fasta', label='Reference genome file')
        known_sites = ListField(
            DataField('variants:vcf', description='Known sites in (.vcf).'),
            label='One or more databases of known polymorphic sites used to exclude regions around known '
                  'polymorphisms from analysis.'
        )
        intervals = DataField(
            data_type='bed',
            label='One or more genomic intervals over which to operate. This field is optional, but it can speed up '
                  'the process by restricting calculations to specific genome regions.'
        )
        read_group = StringField(
            label='Replace read groups in a BAM file.This argument enables the user to replace all read groups in the '
                  'INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM '
                  'file. Addition or replacement is performed using Picard\'s AddOrReplaceReadGroups tool. Input '
                  'should take the form of -name=value delimited by a ``\t``, '
                  'e.g. "-ID=1\t-PL=Illumina\t-SM=sample_1". See [tool\'s documentation]('
                  'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.11.0'
                  '/picard_sam_AddOrReplaceReadGroups.php) for more information on tag names. Note that PL, LB, '
                  'PU and SM are require fields. See caveats of rewriting read groups in the documentation linked '
                  'above.',
            default=''
        )

    class Output:
        """Output fields to BaseQualityScoreRecalibrator."""

        bam = FileField(label='Base quality score recalibrated BAM file')
        bai = FileField(label='Index of base quality score recalibrated BAM file')
        stats = FileField(label='Alignment statistics')
        bigwig = FileField(label='BigWig file')
        species = StringField(label='Species')
        build = StringField(label='Build')
        recal_table = FileField(label='Recalibration tabled')

    def run(self, inputs, outputs):
        """Run the analysis."""
        # Prepare output file names.
        bam = os.path.basename(inputs.bam.bam.path)
        file_name = os.path.splitext(os.path.basename(inputs.bam.bam.path))[0]
        bam_rg = f'{file_name}_RG.bam'

        species = inputs.bam.species

        # Parse read_group argument from a string, delimited by a \t and =
        # into a form that will be accepted by AddOrReplaceReadGroups tool.
        # E.g. '-LB=DAB\t-PL=Illumina\t-PU=barcode\t-SM=sample1' should become
        # ['-LB', 'DAB', '-PL', 'Illumina', '-PU', 'barcode', '-SM', 'sample1']
        # prepended by INPUT and OUTPUT.
        if inputs.read_group:
            arrg = [
                '--INPUT', f'{inputs.bam.bam.path}',
                '--OUTPUT', f'{bam_rg}'
            ]

            for x in inputs.read_group.split('\t'):
                arrg.extend(x.split('='))

            Cmd['gatk']['AddOrReplaceReadGroups'](arrg)
        else:
            shutil.copy2(inputs.bam.bam.path, bam_rg)

        # Make sure the file is indexed.
        Cmd['samtools']['index'](bam_rg)

        recal_table = f'{file_name}_recalibration.table'
        br_inputs = [
            '--input', f'{bam_rg}',
            '--output', f'{recal_table}',
            '--reference', f'{inputs.reference.fasta.path}',
            '--intervals', f'{inputs.intervals.bed.path}'
        ]

        # Add known sites to the input parameters of BaseRecalibrator.
        for site in inputs.known_sites:
            br_inputs.extend(['--known-sites', f'{site.vcf.path}'])

        # Prepare bqsr recalibration file.
        Cmd['gatk']['BaseRecalibrator'](br_inputs)

        self.progress(0.5)

        # Apply base recalibration.
        ab_inputs = [
            '--input', f'{bam_rg}',
            '--reference', f'{inputs.reference.fasta.path}',
            '--bqsr-recal-file', f'{recal_table}',
            '--output', f'{bam}'
        ]
        Cmd['gatk']['ApplyBQSR'](ab_inputs)

        stats = f'{bam}_stats.txt'
        (Cmd['samtools']['flagstat'][f'{bam}'] > stats)()

        self.progress(0.8)

        btb_inputs = [
            f'{bam}',
            f'{species}',
            f'{self.requirements.resources.cores}'
        ]

        Cmd['bamtobigwig.sh'](btb_inputs)

        if not os.path.exists(f'{file_name}.bw'):
            self.info(
                'BigWig file not calculated.'
            )

        self.progress(0.9)

        outputs.bam = bam
        outputs.bai = file_name + '.bai'
        outputs.stats = stats
        outputs.bigwig = file_name + '.bw'
        outputs.species = species
        outputs.build = inputs.bam.build
        outputs.recal_table = recal_table

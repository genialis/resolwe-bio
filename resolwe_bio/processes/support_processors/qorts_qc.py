"""QoRTs QC."""
import os

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    GroupField,
    Process,
    SchedulingClass,
    StringField,
)


class QortsQC(Process):
    """QoRTs QC analysis."""

    slug = 'qorts-qc'
    name = "QoRTs QC"
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/common:1.1.0',
            },
        },
        'resources': {
            'cores': 1,
            'memory': 4096,
        },
    }
    data_name = "QoRTs QC report ({{alignment|sample_name}})"
    version = '1.0.0'
    process_type = 'data:qorts:qc'
    category = 'other'
    entity = {
        'type': 'sample',
        'input': 'alignment',
    }
    scheduling_class = SchedulingClass.BATCH
    description = "Quality of RNA-seq Tool-Set."

    class Input:
        """Input fields."""

        alignment = DataField('alignment:bam', label="Alignment")
        annotation = DataField('annotation:gtf', label="GTF annotation")

        class Options:
            """Options."""

            stranded = StringField(
                label="Assay type",
                default='non_specific',
                choices=[
                    ('non_specific', 'Strand non-specific'),
                    ('forward', 'Strand-specific forward'),
                    ('reverse', 'Strand-specific reverse'),
                ],
            )

        options = GroupField(Options, label="Options")

    class Output:
        """Output fields."""

        plot = FileField(label="QC multiplot")
        summary = FileField(label="QC summary")
        qorts_data = FileField(label="QoRTs report data")

    def run(self, inputs, outputs):
        """Run the analysis."""
        args = [
            'QC',
            '--skipFunctions', 'writeDESeq,writeDEXSeq',
            '--randomSeed', 42,
            '--generatePlots',
            inputs.alignment.bam['file'],
            inputs.annotation.annot['file'],
            'qorts_output',
        ]

        optional_args = []

        # Detect if aligned reads in BAM file are of single or paired-end type
        # The samtools view command counts the number of reads with the SAM flag "read paired (0x1)"
        if Cmd['samtools']('view', '-c', '-f', '1', inputs.alignment.bam['file']).strip() == '0':
            optional_args.append('--singleEnded')

        if inputs.options.stranded == 'forward':
            optional_args.extend(['--stranded', '--stranded_fr_secondstrand'])
        elif inputs.options.stranded == 'reverse':
            optional_args.append('--stranded')

        # Insert optional arguments before the required arguments
        for argument in optional_args:
            args.insert(3, argument)

        Cmd['QoRTs'](args)

        # Compress QoRTs output folder
        Cmd['zip'](['-r', 'qorts_report.zip', 'qorts_output'])

        outputs.plot = 'qorts_output/QC.multiPlot.pdf'
        outputs.summary = 'qorts_output/QC.summary.txt'
        outputs.qorts_data = 'qorts_report.zip'

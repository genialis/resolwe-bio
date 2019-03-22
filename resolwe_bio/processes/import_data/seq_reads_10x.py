"""Import ScRNA-Seq reads."""
import os

from shutil import move
from resolwe.process import (
    FileField,
    ListField,
    Process,
    SchedulingClass,
)


class ImportScRNA10x(Process):
    """Import 10x scRNA reads in FASTQ format."""

    slug = 'upload-sc-10x'
    name = 'Reads (scRNA 10x)'
    process_type = 'data:screads:10x:'
    version = '1.0.0'
    category = 'Import'
    sheduling_class = SchedulingClass.BATCH
    entity = {
        'type': 'sample',
        'descriptor_schema': 'sample',
    }
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/base:ubuntu-18.04'
            }
        }
    }
    data_name = '{{ reads.0.file|default("?") }}'

    class Input:
        """Input fields to process ImportScRNA10x."""

        barcodes = ListField(
            FileField(
                description='Barcodes file(s) in FASTQ format. Usually the forward FASTQ files (R1).',
            ),
            label='Barcodes (.fastq.gz)',
        )
        reads = ListField(
            FileField(
                description='Reads file(s) in FASTQ format. Usually the reverse FASTQ files (R2).',
            ),
            label='Reads (.fastq.gz)',
        )

    class Output:
        """Output fields to process ImportScRNA10x."""

        barcodes = ListField(FileField(), label='Barcodes')
        reads = ListField(FileField(), label='Reads')

    def run(self, inputs, outputs):
        """Run the analysis."""
        # Check if the number of input fastqs is the same
        if len(inputs.barcodes) != len(inputs.reads):
            self.error('The number of reads and barcodes fastqs must be the same.')

        for fastq in inputs.barcodes + inputs.reads:
            move(fastq.file_temp, os.path.basename(fastq.path))

        barcodes_files = [os.path.basename(fastq.path) for fastq in inputs.barcodes]
        reads_files = [os.path.basename(fastq.path) for fastq in inputs.reads]

        outputs.barcodes = barcodes_files
        outputs.reads = reads_files

"""Import ScRNA-Seq reads."""
import os
from shutil import move

from plumbum import TEE

from resolwe.process import (
    Cmd,
    FileField,
    FileHtmlField,
    ListField,
    Process,
    SchedulingClass,
)


class ImportScRNA10x(Process):
    """Import 10x scRNA reads in FASTQ format."""

    slug = "upload-sc-10x"
    name = "Reads (scRNA 10x)"
    process_type = "data:screads:10x:"
    version = "1.2.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.3.1"}},
    }
    data_name = '{{ reads.0.file|default("?") }}'

    class Input:
        """Input fields to process ImportScRNA10x."""

        barcodes = ListField(
            FileField(
                description="Barcodes file(s) in FASTQ format. Usually the forward FASTQ files (R1).",
            ),
            label="Barcodes (.fastq.gz)",
        )
        reads = ListField(
            FileField(
                description="Reads file(s) in FASTQ format. Usually the reverse FASTQ files (R2).",
            ),
            label="Reads (.fastq.gz)",
        )

    class Output:
        """Output fields to process ImportScRNA10x."""

        barcodes = ListField(FileField(), label="Barcodes")
        reads = ListField(FileField(), label="Reads")
        fastqc_url_barcodes = ListField(
            FileHtmlField(), label="Quality control with FastQC (Barcodes)",
        )
        fastqc_url_reads = ListField(
            FileHtmlField(), label="Quality control with FastQC (Reads)",
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        # Check if the number of input fastqs is the same
        if len(inputs.barcodes) != len(inputs.reads):
            self.error("The number of reads and barcodes fastqs must be the same.")

        for fastq in inputs.barcodes + inputs.reads:
            move(fastq.file_temp, os.path.basename(fastq.path))

        barcodes_files = [os.path.basename(fastq.path) for fastq in inputs.barcodes]
        reads_files = [os.path.basename(fastq.path) for fastq in inputs.reads]

        cmd = Cmd["fastqc"]
        for fastq in barcodes_files + reads_files:
            cmd = cmd["{}".format(fastq)]
        cmd = cmd["--extract"]
        cmd = cmd["--outdir=./"]
        _, _, stderr = cmd & TEE
        # FastQC writes both progress and errors to stderr and exits with code 0.
        # Catch if file is empty, wrong format... (Failed to process) or
        # if file path does not exist, file cannot be read (Skipping).
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        barcodes_fastqcs = []
        reads_fastqcs = []
        for barcodes, reads in zip(barcodes_files, reads_files):
            assert barcodes.endswith(".fastq.gz")
            assert reads.endswith(".fastq.gz")
            barcodes_name = barcodes[:-9]
            reads_name = reads[:-9]
            barcodes_fastqcs.append("{}_fastqc.html".format(barcodes_name))
            reads_fastqcs.append("{}_fastqc.html".format(reads_name))

        outputs.barcodes = barcodes_files
        outputs.reads = reads_files
        outputs.fastqc_url_barcodes = barcodes_fastqcs
        outputs.fastqc_url_reads = reads_fastqcs

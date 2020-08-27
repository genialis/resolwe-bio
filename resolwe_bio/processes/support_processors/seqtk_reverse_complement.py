"""Reverse complement reads with Seqtk."""
import os

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    ListField,
    Process,
    StringField,
)


class ReverseComplementSingle(Process):
    """Reverse complement single-end FASTQ reads file using Seqtk."""

    slug = "seqtk-rev-complement-single"
    process_type = "data:reads:fastq:single:seqtk"
    name = "Reverse complement FASTQ (single-end)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "resolwebio/common:1.3.1"},
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    data_name = '{{ reads|sample_name|default("?") }}'
    version = "1.0.1"

    class Input:
        """Input fields to process ReverseComplementSingle."""

        reads = DataField("reads:fastq:single", label="Reads")

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Reverse complemented FASTQ file")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC")
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run the analysis."""
        basename = os.path.basename(inputs.reads.fastq[0].path)
        assert basename.endswith(".fastq.gz")
        name = basename[:-9]
        complemented_name = f"{name}_complemented.fastq"
        # Concatenate multilane reads
        (
            Cmd["cat"][[reads.path for reads in inputs.reads.fastq]]
            > "input_reads.fastq.gz"
        )()

        # Reverse complement reads
        (Cmd["seqtk"]["seq", "-r", "input_reads.fastq.gz"] > complemented_name)()

        _, _, stderr = (
            Cmd["fastqc"][complemented_name, "--extract", "--outdir=./"] & TEE
        )
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        (Cmd["gzip"][complemented_name])()

        outputs.fastq = [f"{complemented_name}.gz"]
        outputs.fastqc_url = [f"{name}_complemented_fastqc.html"]
        outputs.fastqc_archive = [f"{name}_complemented_fastqc.zip"]


class ReverseComplementPaired(Process):
    """Reverse complement paired-end FASTQ reads file using Seqtk."""

    slug = "seqtk-rev-complement-paired"
    process_type = "data:reads:fastq:paired:seqtk"
    name = "Reverse complement FASTQ (paired-end)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "resolwebio/common:1.3.1"},
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    data_name = '{{ reads|sample_name|default("?") }}'
    version = "1.0.1"

    class Input:
        """Input fields to process ReverseComplementPaired."""

        reads = DataField("reads:fastq:paired", label="Reads")
        select_mate = StringField(
            label="Select mate",
            description="Select the which mate should be reverse complemented.",
            choices=[("Mate 1", "Mate 1"), ("Mate 2", "Mate 2"), ("Both", "Both")],
            default="Mate 1",
        )

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Reverse complemented FASTQ file")
        fastq2 = ListField(FileField(), label="Remaining mate")
        fastqc_url = ListField(
            FileHtmlField(), label="Quality control with FastQC (Mate 1)"
        )
        fastqc_archive = ListField(
            FileField(), label="Download FastQC archive (Mate 1)"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Quality control with FastQC (Mate 2)"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download FastQC archive (Mate 2)"
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        basename_mate1 = os.path.basename(inputs.reads.fastq[0].path)
        basename_mate2 = os.path.basename(inputs.reads.fastq2[0].path)
        assert basename_mate1.endswith(".fastq.gz")
        assert basename_mate2.endswith(".fastq.gz")
        name_mate1 = basename_mate1[:-9]
        name_mate2 = basename_mate2[:-9]
        original_mate1 = f"{name_mate1}_original.fastq.gz"
        original_mate2 = f"{name_mate2}_original.fastq.gz"

        (Cmd["cat"][[reads.path for reads in inputs.reads.fastq]] > original_mate1)()
        (Cmd["cat"][[reads.path for reads in inputs.reads.fastq2]] > original_mate2)()

        if inputs.select_mate == "Mate 1":
            complemented_mate1 = f"{name_mate1}_complemented.fastq"
            (Cmd["seqtk"]["seq", "-r", original_mate1] > complemented_mate1)()

            _, _, stderr = (
                Cmd["fastqc"][complemented_mate1, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr or "Skipping" in stderr:
                self.error("Failed while processing with FastQC.")

            _, _, stderr2 = (
                Cmd["fastqc"][original_mate2, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr2 or "Skipping" in stderr2:
                self.error("Failed while processing with FastQC.")

            (Cmd["gzip"][complemented_mate1])()

            outputs.fastq = [f"{complemented_mate1}.gz"]
            outputs.fastq2 = [original_mate2]
            outputs.fastqc_url = [f"{name_mate1}_complemented_fastqc.html"]
            outputs.fastqc_archive = [f"{name_mate1}_complemented_fastqc.zip"]
            outputs.fastqc_url2 = [f"{name_mate2}_original_fastqc.html"]
            outputs.fastqc_archive2 = [f"{name_mate2}_original_fastqc.zip"]

        elif inputs.select_mate == "Mate 2":
            complemented_mate2 = f"{name_mate2}_complemented.fastq"
            (
                Cmd["seqtk"]["seq", "-r", f"{name_mate2}_original.fastq.gz"]
                > complemented_mate2
            )()

            _, _, stderr = (
                Cmd["fastqc"][original_mate1, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr or "Skipping" in stderr:
                self.error("Failed while processing with FastQC.")

            _, _, stderr2 = (
                Cmd["fastqc"][complemented_mate2, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr2 or "Skipping" in stderr2:
                self.error("Failed while processing with FastQC.")

            (Cmd["gzip"][complemented_mate2])()

            outputs.fastq = [original_mate1]
            outputs.fastq2 = [f"{complemented_mate2}.gz"]
            outputs.fastqc_url = [f"{name_mate1}_original_fastqc.html"]
            outputs.fastqc_archive = [f"{name_mate1}_original_fastqc.zip"]
            outputs.fastqc_url2 = [f"{name_mate2}_complemented_fastqc.html"]
            outputs.fastqc_archive2 = [f"{name_mate2}_complemented_fastqc.zip"]

        else:
            complemented_mate1 = f"{name_mate1}_complemented.fastq"
            complemented_mate2 = f"{name_mate2}_complemented.fastq"
            (
                Cmd["seqtk"]["seq", "-r", f"{name_mate1}_original.fastq.gz"]
                > complemented_mate1
            )()

            _, _, stderr = (
                Cmd["fastqc"][complemented_mate1, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr or "Skipping" in stderr:
                self.error("Failed while processing with FastQC.")

            (Cmd["gzip"][complemented_mate1])()

            (
                Cmd["seqtk"]["seq", "-r", f"{name_mate2}_original.fastq.gz"]
                > complemented_mate2
            )()

            _, _, stderr2 = (
                Cmd["fastqc"][complemented_mate2, "--extract", "--outdir=./"] & TEE
            )
            if "Failed to process" in stderr2 or "Skipping" in stderr2:
                self.error("Failed while processing with FastQC.")

            (Cmd["gzip"][complemented_mate2])()

            outputs.fastq = [f"{complemented_mate1}.gz"]
            outputs.fastq2 = [f"{complemented_mate2}.gz"]
            outputs.fastqc_url = [f"{name_mate1}_complemented_fastqc.html"]
            outputs.fastqc_archive = [f"{name_mate1}_complemented_fastqc.zip"]
            outputs.fastqc_url2 = [f"{name_mate2}_complemented_fastqc.html"]
            outputs.fastqc_archive2 = [f"{name_mate2}_complemented_fastqc.zip"]

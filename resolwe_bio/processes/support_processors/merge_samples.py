"""Merge FASTQs."""
from pathlib import Path
from shutil import move

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    ListField,
    Process,
    SchedulingClass,
)


def run_fastqc(fastqs, output_dir):
    """Run fastQC on given FASTQs.

    :param list fastqs: List of fastqs
    :param str output_dir: Output directory

    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    cmd = Cmd["fastqc"]
    for fastq in fastqs:
        cmd = cmd[fastq]
    cmd = cmd["--extract"]
    cmd = cmd[f"--outdir={str(output_path)}"]
    _, _, stderr = cmd & TEE

    return stderr


class MergeFastqSingle(Process):
    """Merge single-end FASTQs into one sample."""

    slug = "merge-fastq-single"
    name = "Merge FASTQ (single-end)"
    process_type = "data:reads:fastq:single"
    version = "1.0.1"
    category = "Other"
    scheduling_class = SchedulingClass.BATCH
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
        "always_create": True,
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.3.1"}},
    }
    data_name = '{{ reads|map("sample_name")|join(", ")|default("?") }}'

    class Input:
        """Input fields to process MergeFastqSingle."""

        reads = ListField(
            DataField(data_type="reads:fastq:single:"), label="Reads data objects",
        )

    class Output:
        """Output fields to process MergeFastqSingle."""

        fastq = ListField(FileField(), label="Reads file")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC",)
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run the analysis."""
        name = f"{Path(inputs.reads[0].fastq[0].path).name.replace('.fastq.gz', '')}_merged"
        merged_fastq = f"{name}.fastq.gz"
        with open(merged_fastq, "wb") as outfile:
            for reads in inputs.reads:
                for fastq in reads.fastq:
                    with open(fastq.path, "rb") as infile:
                        for line in infile:
                            outfile.write(line)

        outputs.fastq = [merged_fastq]

        print("Postprocessing FastQC...")
        stderr = run_fastqc([merged_fastq], "./fastqc")
        # FastQC writes both progress and errors to stderr and exits with code 0.
        # Catch if file is empty, wrong format... (Failed to process) or
        # if file path does not exist, file cannot be read (Skipping).
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        fastqc = f"{name}_fastqc.zip"
        fastqc_path = Path("fastqc")
        move(str(fastqc_path / fastqc), ".")

        report_dir = fastqc_path / f"{name}_fastqc"
        fastqc_url = [
            {
                "file": str(report_dir / "fastqc_report.html"),
                "refs": [str(fastqc_path)],
            }
        ]

        outputs.fastqc_url = fastqc_url
        outputs.fastqc_archive = [fastqc]


class MergeFastqPaired(Process):
    """Merge paired-end FASTQs into one sample."""

    slug = "merge-fastq-paired"
    name = "Merge FASTQ (paired-end)"
    process_type = "data:reads:fastq:paired"
    version = "1.0.1"
    category = "Other"
    scheduling_class = SchedulingClass.BATCH
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
        "always_create": True,
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.3.1"}},
    }
    data_name = '{{ reads|map("sample_name")|join(", ")|default("?") }}'

    class Input:
        """Input fields to process MergeFastqPaired."""

        reads = ListField(
            DataField(data_type="reads:fastq:paired:"), label="Reads data objects",
        )

    class Output:
        """Output fields to process MergeFastqPaired."""

        fastq = ListField(FileField(), label="Reads file (mate 1)")
        fastq2 = ListField(FileField(), label="Reads file (mate 2)")
        fastqc_url = ListField(
            FileHtmlField(), label="Quality control with FastQC (mate 1)",
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Quality control with FastQC (mate 2)",
        )
        fastqc_archive = ListField(
            FileField(), label="Download FastQC archive (mate 1)"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download FastQC archive (mate 2)"
        )

    def run(self, inputs, outputs):
        """Run the analysis."""
        name_1 = f"{Path(inputs.reads[0].fastq[0].path).name.replace('.fastq.gz', '')}_merged"
        name_2 = f"{Path(inputs.reads[0].fastq2[0].path).name.replace('.fastq.gz', '')}_merged"
        merged_fastq_1 = f"{name_1}.fastq.gz"
        merged_fastq_2 = f"{name_2}.fastq.gz"
        with open(merged_fastq_1, "wb") as outfile_1, open(
            merged_fastq_2, "wb"
        ) as outfile_2:
            for reads in inputs.reads:
                for fastq_1, fastq_2 in zip(reads.fastq, reads.fastq2):
                    with open(fastq_1.path, "rb") as infile:
                        for line in infile:
                            outfile_1.write(line)
                    with open(fastq_2.path, "rb") as infile:
                        for line in infile:
                            outfile_2.write(line)

        outputs.fastq = [merged_fastq_1]
        outputs.fastq2 = [merged_fastq_2]

        fastqc_1 = [f"{name_1}_fastqc.zip"]
        fastqc_2 = [f"{name_2}_fastqc.zip"]
        fastqc_path = Path("fastqc")
        report_dir_1 = fastqc_path / f"{name_1}_fastqc"
        fastqc_url_1 = [
            {
                "file": str(report_dir_1 / "fastqc_report.html"),
                "refs": [str(report_dir_1)],
            }
        ]
        report_dir_2 = fastqc_path / f"{name_2}_fastqc"
        fastqc_url_2 = [
            {
                "file": str(report_dir_2 / "fastqc_report.html"),
                "refs": [str(report_dir_2)],
            }
        ]

        print("Postprocessing FastQC...")
        stderr = run_fastqc([merged_fastq_1, merged_fastq_2], "./fastqc")
        # FastQC writes both progress and errors to stderr and exits with code 0.
        # Catch if file is empty, wrong format... (Failed to process) or
        # if file path does not exist, file cannot be read (Skipping).
        if "Failed to process" in stderr or "Skipping" in stderr:
            self.error("Failed while processing with FastQC.")

        for fastqc_zip in Path("fastqc").glob("*_fastqc.zip"):
            move(str(fastqc_zip), ".")

        outputs.fastqc_url = fastqc_url_1
        outputs.fastqc_url2 = fastqc_url_2
        outputs.fastqc_archive = fastqc_1
        outputs.fastqc_archive2 = fastqc_2

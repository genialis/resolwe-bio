"""Run FastQC quality control on paired-end FASTQ reads."""

from pathlib import Path
from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Persistence,
    Process,
    SchedulingClass,
)


class FastQC(Process):
    """Run FastQC on paired-end reads in FASTQ format."""

    slug = "fastqc-paired-end"
    name = "FastQC (paired-end)"
    process_type = "data:fastqc"
    version = "1.0.0"
    category = "QC"
    data_name = '{{ reads.mate1.file|basename|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 2,
            "memory": 4096,
            "storage": 100,
        },
    }

    class Input:
        """Input fields to process FastQC."""

        reads = DataField("reads:fastq:paired", label="Paired-end reads")

    class Output:
        """Output fields to process FastQC."""

        fastqc_report_mate1 = FileField(label="QC report (Mate 1)")
        fastqc_report_mate2 = FileField(label="QC report (Mate 2)")
        fastqc_archive_mate1 = FileField(label="FastQC archive (Mate 1)")
        fastqc_archive_mate2 = FileField(label="FastQC archive (Mate 2)")

    def run(self, inputs, outputs):
        """Run FastQC analysis on paired-end FASTQ files."""

        mates = [
            ("mate1", inputs.reads.output.mate1.path),
            ("mate2", inputs.reads.output.mate2.path),
        ]

        fastqc_reports = {}
        fastqc_archives = {}

        for mate_label, mate_file in mates:
            mate_name = Path(mate_file).name.replace(".fastq.gz", "")
            fastqc_report = mate_name + "_fastqc.html"
            fastqc_archive = mate_name + "_fastqc.zip"
            fastqc_inputs = [
                mate_file,
                "--extract",
                "--outdir=.",
                "-t",
                self.requirements.resources.cores,
            ]
            return_code, _, stderr = Cmd["fastqc"][fastqc_inputs] & TEE(retcode=None)
            if return_code:
                self.error(f"FastQC failed for {mate_label}. ", stderr)

            fastqc_reports[mate_label] = fastqc_report
            fastqc_archives[mate_label] = fastqc_archive

        outputs.fastqc_report_mate1 = fastqc_reports["mate1"]
        outputs.fastqc_report_mate2 = fastqc_reports["mate2"]
        outputs.fastqc_archive_mate1 = fastqc_archives["mate1"]
        outputs.fastqc_archive_mate2 = fastqc_archives["mate2"]

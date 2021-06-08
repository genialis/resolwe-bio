"""Convert aligned reads from BAM to FASTQ format."""
from pathlib import Path

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


class BamToFastqPaired(Process):
    """Convert aligned reads in BAM format to paired-end FASTQ files format."""

    slug = "bamtofastq-paired"
    name = "bedtools bamtofastq (paired-end)"
    category = "BAM processing"
    process_type = "data:reads:fastq:paired:bamtofastq"
    version = "1.0.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.8.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
            "storage": 600,
        },
    }
    entity = {"type": "sample"}
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields for BamToFastqPaired."""

        bam = DataField("alignment:bam", label="BAM file")

    class Output:
        """Output fields for BamToFastqPaired."""

        fastq = ListField(FileField(), label="Remaining mate1 reads")
        fastq2 = ListField(FileField(), label="Remaining mate2 reads")
        fastqc_url = ListField(
            FileHtmlField(), label="Mate1 quality control with FastQC"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Mate2 quality control with FastQC"
        )
        fastqc_archive = ListField(FileField(), label="Download mate1 FastQC archive")
        fastqc_archive2 = ListField(FileField(), label="Download mate2 FastQC archive")

    def run(self, inputs, outputs):
        """Run analysis."""
        name = Path(inputs.bam.output.bam.path).stem
        sorted_bam = f"{name}_sorted.bam"
        mate1 = f"{name}_mate1.fastq"
        mate1_gz = mate1 + ".gz"
        mate2 = f"{name}_mate2.fastq"
        mate2_gz = mate2 + ".gz"

        # For extracted paired-end reads to match, BAM file needs to be name sorted first.
        sort_args = [
            "-@",
            self.requirements.resources.cores,
            "-n",
            "-o",
            sorted_bam,
            inputs.bam.output.bam.path,
        ]
        return_code, _, _ = Cmd["samtools"]["sort"][sort_args] & TEE(retcode=None)
        if return_code:
            self.error("Samtools sort command failed.")

        self.progress(0.3)

        # Convert aligned reads into paired-end FASTQ files
        extract_args = ["-i", sorted_bam, "-fq", mate1, "-fq2", mate2]
        (Cmd["bedtools"]["bamtofastq"][extract_args] > "log.txt")()

        self.progress(0.8)

        (Cmd["gzip"][mate1])()
        (Cmd["gzip"][mate2])()

        self.progress(0.9)

        # Prepare final FASTQC report
        fastqc_args_mate1 = [
            mate1_gz,
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]
        return_code, _, _ = Cmd["fastqc.sh"][fastqc_args_mate1] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing FASTQC report.")

        self.progress(0.95)

        fastqc_args_mate2 = [
            mate2_gz,
            "fastqc",
            "fastqc_archive2",
            "fastqc_url2",
        ]
        return_code, _, _ = Cmd["fastqc.sh"][fastqc_args_mate2] & TEE(retcode=None)
        if return_code:
            self.error("Error while preparing FASTQC report.")

        # Save the outputs
        outputs.fastq = [mate1_gz]
        outputs.fastq2 = [mate2_gz]

"""UMI-tools dedup."""

import os
from glob import glob

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Process,
    SchedulingClass,
    StringField,
)


class UmiToolsDedup(Process):
    """Deduplicate reads using UMI and mapping coordinates."""

    slug = "umi-tools-dedup"
    name = "UMI-tools dedup"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0",
            },
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    data_name = "{{ alignment|name|default('?') }}"
    version = "1.5.1"
    process_type = "data:alignment:bam:umitools:dedup"
    category = "FASTQ processing"
    entity = {
        "type": "sample",
        "input": "alignment",
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        alignment = DataField("alignment:bam", label="Alignment")

    class Output:
        """Output fields."""

        bam = FileField(label="Clipped BAM file")
        bai = FileField(label="Index of clipped BAM file")
        stats = FileField(label="Alignment statistics")
        dedup_log = FileField(label="Deduplication log")
        dedup_stats = FileField(label="Deduplication stats")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        alignment_path = os.path.basename(inputs.alignment.output.bam.path)
        assert alignment_path.endswith(".bam")
        name = alignment_path[:-4]

        out_bam = "{}_dedup.bam".format(name)
        out_log = "{}_dedup.log".format(name)
        stats = "{}_stats.txt".format(name)
        dedup_stats = "{}_dedup_stats.zip".format(name)

        args = [
            "dedup",
            "-I",
            inputs.alignment.output.bam.path,
            "-S",
            out_bam,
            "-L",
            out_log,
            "--multimapping-detection-method=NH",
            "--output-stats=dedup_stats",
        ]

        # Detect if aligned reads in BAM file are of single or paired-end type
        # The samtools view command counts the number of reads with the SAM flag "read paired (0x1)"
        if (
            Cmd["samtools"](
                "view", "-c", "-f", "1", inputs.alignment.output.bam.path
            ).strip()
            != "0"
        ):
            args.append("--paired")

        # Run UMI-tools dedup
        return_code, _, _ = Cmd["umi_tools"][args] & TEE(retcode=None)
        if return_code:
            self.error("Deduplication of {}.bam failed.".format(name))

        # Compress deduplication stats files
        Cmd["zip"]([dedup_stats, *glob("dedup_stats_*")])

        # Index deduplicated output .bam file
        return_code, _, _ = Cmd["samtools"]["index", out_bam] & TEE(retcode=None)
        if return_code:
            self.error("Indexing of {} failed.".format(out_bam))

        # Calculate alignment statistics
        (Cmd["samtools"]["flagstat", out_bam] > stats)()

        # Save the outputs
        outputs.bam = out_bam
        outputs.bai = out_bam + ".bai"
        outputs.stats = stats
        outputs.dedup_log = out_log
        outputs.dedup_stats = dedup_stats
        outputs.species = inputs.alignment.output.species
        outputs.build = inputs.alignment.output.build

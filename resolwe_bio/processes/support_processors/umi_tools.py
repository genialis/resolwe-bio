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
                "image": "resolwebio/rnaseq:4.9.0",
            },
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    data_name = "UMI-tools dedup ({{alignment|sample_name}})"
    version = "1.1.1"
    process_type = "data:alignment:bam:umitools:dedup"
    category = "Other"
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
        bigwig = FileField(label="BigWig file", required=False)
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        alignment_path = os.path.basename(inputs.alignment.bam.path)
        assert alignment_path.endswith(".bam")
        name = alignment_path[:-4]

        out_bam = "{}_dedup.bam".format(name)
        out_log = "{}_dedup.log".format(name)
        stats = "{}_stats.txt".format(name)
        dedup_stats = "{}_dedup_stats.zip".format(name)

        args = [
            "dedup",
            "-I",
            inputs.alignment.bam.path,
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
            Cmd["samtools"]("view", "-c", "-f", "1", inputs.alignment.bam.path).strip()
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

        # Calculate BigWig file
        bigwig_args = [
            out_bam,
            inputs.alignment.species,
            self.requirements.resources.cores,
        ]
        return_code, _, _ = Cmd["bamtobigwig.sh"][bigwig_args] & TEE(retcode=None)
        if return_code:
            self.error("Failed to calculate BigWig file.")

        # Save the outputs
        outputs.bam = out_bam
        outputs.bai = out_bam + ".bai"
        outputs.stats = stats
        outputs.dedup_log = out_log
        outputs.dedup_stats = dedup_stats
        outputs.species = inputs.alignment.species
        outputs.build = inputs.alignment.build

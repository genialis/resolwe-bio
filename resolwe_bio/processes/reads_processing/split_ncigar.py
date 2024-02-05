"""Run GATK SplitNCigarReads."""

import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class GatkSplitNCigarReads(Process):
    """Splits reads that contain Ns in their cigar string.

    Identifies all N cigar elements and creates k+1 new reads (where k is the number
    of N cigar elements). The first read includes the bases that are to the left of
    the first N element, while the part of the read that is to the right of the N
    (including the Ns) is hard clipped and so on for the rest of the new reads. Used
    for post-processing RNA reads aligned against the full reference.
    """

    slug = "gatk-split-ncigar"
    name = "GATK SplitNCigarReads"
    category = "GATK"
    process_type = "data:alignment:bam:splitncigar"
    version = "1.2.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
            "storage": 200,
        },
    }
    entity = {"type": "sample"}

    data_name = "{{ bam|name|default('?') }}"

    class Input:
        """Input fields for GatkSplitNCigarReads."""

        bam = DataField(
            data_type="alignment:bam",
            label="Alignment BAM file",
        )
        ref_seq = DataField(
            data_type="seq:nucleotide",
            label="Reference sequence FASTA file",
        )

        class Advanced:
            """Advanced options."""

            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields for GatkSplitNCigarReads."""

        bam = FileField(label="BAM file with reads split at N CIGAR elements")
        bai = FileField(label="Index of BAM file")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        file_name = Path(inputs.bam.output.bam.path).name[:-4]
        bam = file_name + ".splitNcigar.bam"
        bai = file_name + ".splitNcigar.bai"

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-I",
            inputs.bam.output.bam.path,
            "-O",
            bam,
            "--tmp-dir",
            TMPDIR,
        ]

        return_code, stdout, stderr = Cmd["gatk"]["SplitNCigarReads"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error(
                "GATK SplitNCigarReads failed. Check standard output for more details."
            )

        stats = f"{bam}_stats.txt"
        (Cmd["samtools"]["flagstat"][f"{bam}"] > stats)()

        outputs.bam = bam
        outputs.bai = bai
        outputs.stats = stats
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build

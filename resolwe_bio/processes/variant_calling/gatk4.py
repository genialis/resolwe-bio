"""Run GATK HaplotypeCaller in GVCF mode."""

import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    Process,
    SchedulingClass,
    StringField,
)


class GatkHaplotypeCallerGvcf(Process):
    """Run GATK HaplotypeCaller in GVCF mode."""

    slug = "gatk-haplotypecaller-gvcf"
    name = "GATK HaplotypeCaller (GVCF)"
    category = "GATK"
    process_type = "data:variants:gvcf"
    version = "1.3.0"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
        },
    }
    data_name = "{{ bam|name|default('?') }}"

    class Input:
        """Input fields for GatkHaplotypeCallerGvcf."""

        bam = DataField("alignment:bam", label="Analysis ready BAM file")
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

        class Options:
            """Options."""

            intervals = DataField(
                "bed",
                label="Use intervals BED file to limit the analysis to the specified parts of the genome.",
                required=False,
            )

            contamination = FloatField(
                label="Contamination fraction",
                default=0,
                description="Fraction of contamination in sequencing data (for all samples) to aggressively remove.",
            )

        options = GroupField(Options, label="Options")

    class Output:
        """Output fields for GatkHaplotypeCallerGvcf."""

        vcf = FileField(label="GVCF file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        name = Path(inputs.bam.output.bam.path).stem
        variants = name + ".g.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"

        args = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-I",
            inputs.bam.output.bam.path,
            "-O",
            variants,
            "-contamination",
            inputs.options.contamination,
            "--tmp-dir",
            TMPDIR,
            "-G",
            "StandardAnnotation",
            "-G",
            "StandardHCAnnotation",
            "-G",
            "AS_StandardAnnotation",
            "-GQB",
            10,
            "-GQB",
            20,
            "-GQB",
            30,
            "-GQB",
            40,
            "-GQB",
            50,
            "-GQB",
            60,
            "-GQB",
            70,
            "-GQB",
            80,
            "-GQB",
            90,
            "-ERC",
            "GVCF",
        ]

        if inputs.options.intervals:
            args.extend(["-L", inputs.options.intervals.output.bed.path])

        return_code, _, _ = Cmd["gatk"]["HaplotypeCaller"][args] & TEE(retcode=None)
        if return_code:
            self.error("GATK HaplotypeCaller tool failed.")

        # Compress and index the output variants file
        (Cmd["bgzip"]["-c", variants] > variants_gz)()
        Cmd["tabix"]["-p", "vcf", variants_gz]()

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build

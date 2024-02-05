"""Run GATK HaplotypeCaller."""

import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class GatkHaplotypeCaller(Process):
    """GATK HaplotypeCaller Variant Calling.

    Call germline SNPs and indels via local re-assembly of haplotypes.

    The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local
    de-novo assembly of haplotypes in an active region. In other words, whenever the program
    encounters a region showing signs of variation, it discards the existing mapping information
    and completely reassembles the reads in that region. This allows the HaplotypeCaller to be
    more accurate when calling regions that are traditionally difficult to call, for example when
    they contain different types of variants close to each other. It also makes the HaplotypeCaller
    much better at calling indels than position-based callers like UnifiedGenotyper.
    """

    slug = "vc-gatk4-hc"
    name = "GATK4 (HaplotypeCaller)"
    category = "GATK"
    process_type = "data:variants:vcf:gatk:hc"
    version = "1.5.0"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    data_name = "{{ alignment|name|default('?') }}"

    class Input:
        """Input fields for GatkHaplotypeCaller."""

        alignment = DataField(
            data_type="alignment:bam", label="Analysis ready BAM file"
        )
        genome = DataField(data_type="seq:nucleotide", label="Reference genome")
        intervals_bed = DataField(
            data_type="bed",
            label="Intervals (from BED file)",
            description="Use this option to perform the analysis over only part of the genome.",
            required=False,
        )
        dbsnp = DataField(
            data_type="variants:vcf",
            label="dbSNP file",
            description="Database of known polymorphic sites.",
        )
        stand_call_conf = IntegerField(
            label="Min call confidence threshold",
            default=30,
            description="The minimum phred-scaled confidence threshold at which "
            "variants should be called.",
        )
        mbq = IntegerField(
            label="Min Base Quality",
            default=20,
            description="Minimum base quality required to consider a base for calling.",
        )
        max_reads = IntegerField(
            label="Max reads per aligment start site",
            default=50,
            description="Maximum number of reads to retain per alignment start position. "
            "Reads above this threshold will be downsampled. Set to 0 to disable.",
        )

        class Advanced:
            """Advanced options."""

            interval_padding = IntegerField(
                label="Interval padding",
                required=False,
                description="Amount of padding (in bp) to add to each interval "
                "you are including. The recommended value is 100.",
                hidden="!intervals_bed",
            )
            soft_clipped = BooleanField(
                label="Do not analyze soft clipped bases in the reads",
                default=False,
                description="Suitable option for RNA-seq variant calling.",
            )
            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of "
                "the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields for GatkHaplotypeCaller."""

        vcf = FileField(label="VCF file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        name = Path(inputs.alignment.output.bam.path).stem
        variants = name + ".gatkHC.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "-R",
            inputs.genome.output.fasta.path,
            "-I",
            inputs.alignment.output.bam.path,
            "-O",
            variants,
            "--dbsnp",
            inputs.dbsnp.output.vcf.path,
            "--min-base-quality-score",
            inputs.mbq,
            "--max-reads-per-alignment-start",
            inputs.max_reads,
            "--standard-min-confidence-threshold-for-calling",
            inputs.stand_call_conf,
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.advanced.soft_clipped:
            args.append("--dont-use-soft-clipped-bases")

        if inputs.intervals_bed:
            args.extend(["-L", inputs.intervals_bed.output.bed.path])

            if inputs.advanced.interval_padding:
                args.extend(["--interval-padding", inputs.advanced.interval_padding])

        return_code, stdout, stderr = Cmd["gatk"]["HaplotypeCaller"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("GATK HaplotypeCaller tool failed.")

        self.progress(0.8)

        # Compress and index the output variants file
        (Cmd["bgzip"]["-c", variants] > variants_gz)()
        self.progress(0.9)

        Cmd["tabix"]["-p", "vcf", variants_gz]()
        self.progress(0.95)

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = inputs.alignment.output.species
        outputs.build = inputs.alignment.output.build

"""Run GATK SelectVariants."""

import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


def return_sample_count(vcf, error):
    """Count number of samples in the input VCF file."""
    try:
        return int((Cmd["bcftools"]["query"]["-l", vcf] | Cmd["wc"]["-l"])().strip())
    except Exception as err:
        error(
            f"Unable to determine sample count in VCF file. Original error was: {err}"
        )


class GatkSelectVariants(Process):
    """Select a subset of variants based on various criteria using GATK SelectVariants.

    This tool works with multi-sample VCF file as an input.
    """

    slug = "gatk-select-variants"
    process_type = "data:variants:vcf:selectvariants"
    name = "GATK SelectVariants (multi-sample)"
    version = "1.2.0"
    category = "GATK"
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

    data_name = "Selected variants"

    class Input:
        """Input fields for GatkSelectVariants."""

        vcf = DataField("variants:vcf", label="Input data (VCF)")

        intervals = DataField(
            "bed",
            label="Intervals file (.bed)",
            description="One or more genomic intervals over which to operate. This can also be "
            "used to get data from a specific interval.",
            required=False,
        )

        select_type = ListField(
            StringField(),
            label="Select only a certain type of variants from the input file",
            description="This argument selects particular kinds of variants out of a list. If "
            "left empty, there is no type selection and all variant types are considered for "
            "other selection criteria. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, "
            "NO_VARIATION. Can be specified multiple times.",
            required=False,
        )

        exclude_filtered = BooleanField(
            label="Don't include filtered sites",
            default=False,
            description="If this flag is enabled, sites that have been marked as filtered (i.e. have "
            "anything other than `.` or `PASS` in the FILTER field) will be excluded from the output.",
        )

        class AdvancedOptions:
            """Advanced options."""

            ref_seq = DataField(
                "seq:nucleotide",
                label="Reference sequence",
                required=False,
            )

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

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GatkSelectVariants."""

        vcf = FileField(label="Selected variants (VCF)")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        # check the VCF file content
        sample_count = return_sample_count(
            vcf=inputs.vcf.output.vcf.path, error=self.error
        )
        if not sample_count > 1:
            self.error(
                f"The input VCF file should contain data for multiple samples. "
                f"The input contains data for {sample_count} sample(s)."
            )

        selected_variants = "selected_variants.vcf.gz"
        selected_variants_index = selected_variants + ".tbi"
        species = inputs.vcf.output.species
        build = inputs.vcf.output.build

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced_options.java_gc_threads
        )
        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced_options.max_heap_size}g",
            "--variant",
            inputs.vcf.output.vcf.path,
            "--output",
            selected_variants,
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.advanced_options.ref_seq:
            if (
                inputs.advanced_options.ref_seq.output.species
                != inputs.vcf.output.species
            ):
                self.error(
                    "The species information of the provided reference "
                    "sequence file does not match the species of the input VCF."
                )
            if inputs.advanced_options.ref_seq.output.build != inputs.vcf.output.build:
                self.error(
                    "The genome build information of the provided reference "
                    "sequence file does not match the build of the input VCFs."
                )
            args.extend(["-R", inputs.advanced_options.ref_seq.output.fasta.path])

        if inputs.intervals:
            args.extend(["-L", inputs.intervals.output.bed.path])
        if inputs.select_type:
            for type in inputs.select_type:
                args.extend(["-select-type", type])

        if inputs.exclude_filtered:
            args.append("--exclude-filtered")

        return_code, stdout, stderr = Cmd["gatk"]["SelectVariants"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("GATK SelectVariants failed.")

        outputs.vcf = selected_variants
        outputs.tbi = selected_variants_index
        outputs.species = species
        outputs.build = build


class GatkSelectVariantsSingleSample(Process):
    """Select a subset of variants based on various criteria using GATK SelectVariants.

    This tool works with single-sample VCF file as an input.
    """

    slug = "gatk-select-variants-single"
    process_type = "data:variants:vcf:selectvariants:single"
    name = "GATK SelectVariants (single-sample)"
    version = "1.1.0"
    entity = {
        "type": "sample",
    }
    category = "GATK"
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

    data_name = "{{ vcf|name|default('?') }}"

    class Input:
        """Input fields for GatkSelectVariantsSingleSample."""

        vcf = DataField("variants:vcf", label="Input data (VCF)")

        intervals = DataField(
            "bed",
            label="Intervals file (.bed)",
            description="One or more genomic intervals over which to operate. This can also be "
            "used to get data from a specific interval.",
            required=False,
        )

        select_type = ListField(
            StringField(),
            label="Select only a certain type of variants from the input file",
            description="This argument selects particular kinds of variants out of a list. If "
            "left empty, there is no type selection and all variant types are considered for "
            "other selection criteria. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, "
            "NO_VARIATION. Can be specified multiple times.",
            required=False,
        )

        exclude_filtered = BooleanField(
            label="Don't include filtered sites",
            default=False,
            description="If this flag is enabled, sites that have been marked as filtered (i.e. have "
            "anything other than `.` or `PASS` in the FILTER field) will be excluded from the output.",
        )

        class AdvancedOptions:
            """Advanced options."""

            ref_seq = DataField(
                "seq:nucleotide",
                label="Reference sequence",
                required=False,
            )

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

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GatkSelectVariantsSingleSample."""

        vcf = FileField(label="Selected variants (VCF)")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        # check the VCF file content
        sample_count = return_sample_count(
            vcf=inputs.vcf.output.vcf.path, error=self.error
        )
        if sample_count != 1:
            self.error(
                f"The input VCF should contain data for a single sample. "
                f"The input contains data for {sample_count} sample(s)."
            )

        selected_variants = "selected_variants.vcf.gz"
        selected_variants_index = selected_variants + ".tbi"
        species = inputs.vcf.output.species
        build = inputs.vcf.output.build

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced_options.java_gc_threads
        )
        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced_options.max_heap_size}g",
            "--variant",
            inputs.vcf.output.vcf.path,
            "--output",
            selected_variants,
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.advanced_options.ref_seq:
            if (
                inputs.advanced_options.ref_seq.output.species
                != inputs.vcf.output.species
            ):
                self.error(
                    "The species information of the provided reference "
                    "sequence file does not match the species of the input VCF."
                )
            if inputs.advanced_options.ref_seq.output.build != inputs.vcf.output.build:
                self.error(
                    "The genome build information of the provided reference "
                    "sequence file does not match the build of the input VCF."
                )
            args.extend(["-R", inputs.advanced_options.ref_seq.output.fasta.path])

        if inputs.intervals:
            args.extend(["-L", inputs.intervals.output.bed.path])
        if inputs.select_type:
            for type in inputs.select_type:
                args.extend(["-select-type", type])

        if inputs.exclude_filtered:
            args.append("--exclude-filtered")

        return_code, stdout, stderr = Cmd["gatk"]["SelectVariants"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("GATK SelectVariants failed.")

        outputs.vcf = selected_variants
        outputs.tbi = selected_variants_index
        outputs.species = species
        outputs.build = build

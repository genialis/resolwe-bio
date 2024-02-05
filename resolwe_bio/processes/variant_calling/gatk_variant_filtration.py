"""Run GATK VariantFiltration."""

import os

from plumbum import TEE

from resolwe.process import (
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


class GatkVariantFiltration(Process):
    """Filter multi-sample variant calls based on INFO and/or FORMAT annotations.

    This tool is designed for hard-filtering variant calls based on certain criteria.
    Records are hard-filtered by changing the value in the FILTER field to something
    other than PASS. Passing variants are annotated as PASS and failing variants are
    annotated with the name(s) of the filter(s) they failed. If you want to remove
    failing variants, use GATK SelectVariants process.
    """

    slug = "gatk-variant-filtration"
    process_type = "data:variants:vcf:variantfiltration"
    name = "GATK VariantFiltration (multi-sample)"
    version = "1.3.0"
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
    data_name = "Filtered variants"

    class Input:
        """Input fields for GatkVariantFiltration."""

        vcf = DataField(data_type="variants:vcf", label="Input data (VCF)")
        ref_seq = DataField(data_type="seq:nucleotide", label="Reference sequence")
        filter_expressions = ListField(
            StringField(),
            label="Expressions used with INFO fields to filter",
            description="VariantFiltration accepts any number of JEXL expressions "
            "(so you can have two named filters by using --filter-name One "
            "--filter-expression 'X < 1' --filter-name Two --filter-expression 'X > 2'). "
            "It is preferable to use multiple expressions, each specifying an individual "
            "filter criteria, to a single compound expression that specifies multiple "
            "filter criteria. Input expressions one by one and press ENTER after each "
            "expression. Examples of filter expression: 'FS > 30', 'DP > 10'.",
            required=False,
        )
        filter_name = ListField(
            StringField(),
            label="Names to use for the list of filters",
            description="This name is put in the FILTER field for variants that get "
            "filtered. Note that there must be a 1-to-1 mapping between filter expressions "
            "and filter names. Input expressions one by one and press ENTER after each name. "
            "Warning: filter names should be in the same order as filter expressions. "
            "Example: you specified filter expressions 'FS > 30' and 'DP > 10', now "
            "specify filter names 'FS' and 'DP'.",
            required=False,
        )
        genotype_filter_expressions = ListField(
            StringField(),
            label="Expressions used with FORMAT field to filter",
            description="Similar to the INFO field based expressions, but used on the FORMAT "
            "(genotype) fields instead. VariantFiltration will add the sample-level FT tag to "
            "the FORMAT field of filtered samples (this does not affect the record's FILTER tag). "
            "One can filter normally based on most fields (e.g. 'GQ < 5.0'), but the GT "
            "(genotype) field is an exception. We have put in convenience methods so that "
            "one can now filter out hets ('isHet == 1'), refs ('isHomRef == 1'), or homs "
            "('isHomVar == 1'). Also available are expressions isCalled, isNoCall, isMixed, "
            "and isAvailable, in accordance with the methods of the Genotype object. "
            "To filter by alternative allele depth, use the expression: 'AD.1 < 5'. This "
            "filter expression will filter all the samples in the multi-sample VCF file.",
            required=False,
        )
        genotype_filter_name = ListField(
            StringField(),
            label="Names to use for the list of genotype filters",
            description="Similar to the INFO field based expressions, but used on the FORMAT "
            "(genotype) fields instead. Warning: filter names should be in the same order as "
            "filter expressions.",
            required=False,
        )
        mask = DataField(
            data_type="variants:vcf",
            label="Input mask",
            description="Any variant which overlaps entries from the provided "
            "mask file will be filtered.",
            required=False,
        )
        mask_name = StringField(
            label="The text to put in the FILTER field if a 'mask' is provided",
            description="When using the mask file, the mask name will be annotated in "
            "the variant record.",
            required=False,
            disabled="!mask",
        )

        class Advanced:
            """Advanced options."""

            cluster = IntegerField(
                label="Cluster size",
                default=3,
                description="The number of SNPs which make up a cluster. Must be at least 2.",
            )
            window = IntegerField(
                label="Window size",
                default=0,
                description="The window size (in bases) in which to evaluate clustered SNPs.",
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

        advanced = GroupField(
            Advanced,
            label="Advanced options",
        )

    class Output:
        """Output fields for GatkVariantFiltration."""

        vcf = FileField(label="Filtered variants (VCF)")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        filtered_variants = "filtered_variants.vcf.gz"
        filtered_variants_index = filtered_variants + ".tbi"

        # check the VCF file content
        sample_count = return_sample_count(
            vcf=inputs.vcf.output.vcf.path, error=self.error
        )
        if not sample_count > 1:
            self.error(
                f"The input VCF file should contain data for multiple samples. "
                f"The input contains data for {sample_count} sample(s)."
            )

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "-V",
            inputs.vcf.output.vcf.path,
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-O",
            filtered_variants,
            "--window",
            inputs.advanced.window,
            "--cluster",
            inputs.advanced.cluster,
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.filter_expressions:
            if len(inputs.filter_expressions) != len(inputs.filter_name):
                self.error(
                    "The number of filter expressions and filter names is not the same."
                )
            for name, exp in zip(inputs.filter_name, inputs.filter_expressions):
                args.extend(["--filter-name", name, "--filter-expression", exp])

        if inputs.genotype_filter_expressions:
            if len(inputs.genotype_filter_expressions) != len(
                inputs.genotype_filter_name
            ):
                self.error(
                    "The number of genotype filter expressions and filter names is not the same."
                )

            for name, exp in zip(
                inputs.genotype_filter_name, inputs.genotype_filter_expressions
            ):
                args.extend(
                    [
                        "--genotype-filter-name",
                        name,
                        "--genotype-filter-expression",
                        exp,
                    ]
                )

        if inputs.mask:
            if not inputs.mask_name:
                self.error(
                    "If you specify a mask file, please specify 'mask name' - the text to "
                    "put in the FILTER field"
                )
            args.extend(
                ["--mask", inputs.mask.output.vcf.path, "--mask-name", inputs.mask_name]
            )

        return_code, stdout, stderr = Cmd["gatk"]["VariantFiltration"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error(
                "GATK VariantFiltration failed. Check standard output for more "
                "information."
            )

        outputs.vcf = filtered_variants
        outputs.tbi = filtered_variants_index
        outputs.species = inputs.vcf.output.species
        outputs.build = inputs.vcf.output.build


class GatkVariantFiltrationSingle(Process):
    """Filter single-sample variant calls based on INFO and/or FORMAT annotations.

    This tool is designed for hard-filtering variant calls based on certain criteria.
    Records are hard-filtered by changing the value in the FILTER field to something
    other than PASS. Passing variants are annotated as PASS and failing variants are
    annotated with the name(s) of the filter(s) they failed. If you want to remove
    failing variants, use GATK SelectVariants process.
    """

    slug = "gatk-variant-filtration-single"
    process_type = "data:variants:vcf:variantfiltration:single"
    name = "GATK VariantFiltration (single-sample)"
    version = "1.3.0"
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
        """Input fields for GatkVariantFiltrationSingle."""

        vcf = DataField(data_type="variants:vcf", label="Input data (VCF)")
        ref_seq = DataField(data_type="seq:nucleotide", label="Reference sequence")
        filter_expressions = ListField(
            StringField(),
            label="Expressions used with INFO fields to filter",
            description="VariantFiltration accepts any number of JEXL expressions "
            "(so you can have two named filters by using --filter-name One "
            "--filter-expression 'X < 1' --filter-name Two --filter-expression 'X > 2'). "
            "It is preferable to use multiple expressions, each specifying an individual "
            "filter criteria, to a single compound expression that specifies multiple "
            "filter criteria. Input expressions one by one and press ENTER after each "
            "expression. Examples of filter expression: 'FS > 30', 'DP > 10'.",
            required=False,
        )
        filter_name = ListField(
            StringField(),
            label="Names to use for the list of filters",
            description="This name is put in the FILTER field for variants that get "
            "filtered. Note that there must be a 1-to-1 mapping between filter expressions "
            "and filter names. Input expressions one by one and press ENTER after each name. "
            "Warning: filter names should be in the same order as filter expressions. "
            "Example: you specified filter expressions 'FS > 30' and 'DP > 10', now "
            "specify filter names 'FS' and 'DP'.",
            required=False,
        )
        genotype_filter_expressions = ListField(
            StringField(),
            label="Expressions used with FORMAT field to filter",
            description="Similar to the INFO field based expressions, but used on the FORMAT "
            "(genotype) fields instead. VariantFiltration will add the sample-level FT tag to "
            "the FORMAT field of filtered samples (this does not affect the record's FILTER tag). "
            "One can filter normally based on most fields (e.g. 'GQ < 5.0'), but the GT "
            "(genotype) field is an exception. We have put in convenience methods so that "
            "one can now filter out hets ('isHet == 1'), refs ('isHomRef == 1'), or homs "
            "('isHomVar == 1'). Also available are expressions isCalled, isNoCall, isMixed, "
            "and isAvailable, in accordance with the methods of the Genotype object. "
            "To filter by alternative allele depth, use the expression: 'AD.1 < 5'.",
            required=False,
        )
        genotype_filter_name = ListField(
            StringField(),
            label="Names to use for the list of genotype filters",
            description="Similar to the INFO field based expressions, but used on the FORMAT "
            "(genotype) fields instead. Warning: filter names should be in the same order as "
            "filter expressions.",
            required=False,
        )
        mask = DataField(
            data_type="variants:vcf",
            label="Input mask",
            description="Any variant which overlaps entries from the provided "
            "mask file will be filtered.",
            required=False,
        )
        mask_name = StringField(
            label="The text to put in the FILTER field if a 'mask' is provided",
            description="When using the mask file, the mask name will be annotated in "
            "the variant record.",
            required=False,
            disabled="!mask",
        )

        class Advanced:
            """Advanced options."""

            cluster = IntegerField(
                label="Cluster size",
                default=3,
                description="The number of SNPs which make up a cluster. Must be at least 2.",
            )
            window = IntegerField(
                label="Window size",
                default=0,
                description="The window size (in bases) in which to evaluate clustered SNPs.",
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

        advanced = GroupField(
            Advanced,
            label="Advanced options",
        )

    class Output:
        """Output fields for GatkVariantFiltrationSingle."""

        vcf = FileField(label="Filtered variants (VCF)")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        filtered_variants = "filtered_variants.vcf.gz"
        filtered_variants_index = filtered_variants + ".tbi"

        # check the VCF file content
        sample_count = return_sample_count(
            vcf=inputs.vcf.output.vcf.path, error=self.error
        )
        if sample_count != 1:
            self.error(
                f"The input VCF should contain data for a single sample. "
                f"The input contains data for {sample_count} sample(s)."
            )

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        args = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "-V",
            inputs.vcf.output.vcf.path,
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-O",
            filtered_variants,
            "--window",
            inputs.advanced.window,
            "--cluster",
            inputs.advanced.cluster,
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.filter_expressions:
            if len(inputs.filter_expressions) != len(inputs.filter_name):
                self.error(
                    "The number of filter expressions and filter names is not the same."
                )
            for name, exp in zip(inputs.filter_name, inputs.filter_expressions):
                args.extend(["--filter-name", name, "--filter-expression", exp])

        if inputs.genotype_filter_expressions:
            if len(inputs.genotype_filter_expressions) != len(
                inputs.genotype_filter_name
            ):
                self.error(
                    "The number of genotype filter expressions and filter names is not the same."
                )

            for name, exp in zip(
                inputs.genotype_filter_name, inputs.genotype_filter_expressions
            ):
                args.extend(
                    [
                        "--genotype-filter-name",
                        name,
                        "--genotype-filter-expression",
                        exp,
                    ]
                )

        if inputs.mask:
            if not inputs.mask_name:
                self.error(
                    "If you specify a mask file, please specify 'mask name' - the text to "
                    "put in the FILTER field"
                )
            args.extend(
                ["--mask", inputs.mask.output.vcf.path, "--mask-name", inputs.mask_name]
            )

        return_code, stdout, stderr = Cmd["gatk"]["VariantFiltration"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error(
                "GATK VariantFiltration failed. Check standard output for more "
                "information."
            )

        outputs.vcf = filtered_variants
        outputs.tbi = filtered_variants_index
        outputs.species = inputs.vcf.output.species
        outputs.build = inputs.vcf.output.build

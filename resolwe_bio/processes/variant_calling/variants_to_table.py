"""Run GATK VariantsToTable."""

import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


class GatkVariantsToTable(Process):
    """Run GATK VariantsToTable.

    This tool extracts specified fields for each variant in a VCF file
    to a tab-delimited table, which may be easier to work with than a VCF.
    For additional information, please see
    [manual page](https://gatk.broadinstitute.org/hc/en-us/articles/360036711531-VariantsToTable)
    """

    slug = "variants-to-table"
    name = "GATK VariantsToTable"
    category = "GATK"
    process_type = "data:variantstable"
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
    data_name = "Variants in table"

    class Input:
        """Input fields for GATK VariantsToTable."""

        vcf = DataField("variants:vcf", label="Input VCF file")

        vcf_fields = ListField(
            StringField(),
            label="Select VCF fields",
            description="The name of a standard VCF field or an "
            "INFO field to include in the output table. "
            "The field can be any standard VCF column (e.g. CHROM, ID, QUAL) "
            "or any annotation name in the INFO field (e.g. AC, AF).",
            default=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
            ],
        )

        class AdvancedOptions:
            """Advanced options."""

            gf_fields = ListField(
                StringField(),
                label="Include FORMAT/sample-level fields",
                default=[
                    "GT",
                    "GQ",
                ],
            )
            split_alleles = BooleanField(
                label="Split multi-allelic records into multiple lines",
                description="By default, a variant record with multiple "
                "ALT alleles will be summarized in one line, with per "
                "alt-allele fields (e.g. allele depth) separated by commas."
                "This may cause difficulty when the table is loaded by "
                "an R script, for example. Use this flag to write multi-allelic "
                "records on separate lines of output.",
                default=True,
            )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for GATK VariantsToTable."""

        tsv = FileField(label="Tab-delimited file with variants")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        variants_table = "variants_table.tsv"

        args = [
            "-V",
            inputs.vcf.output.vcf.path,
            "-O",
            variants_table,
            "--tmp-dir",
            TMPDIR,
        ]

        for field in inputs.vcf_fields:
            args.extend(["-F", field])

        for gf_field in inputs.advanced_options.gf_fields:
            args.extend(["-GF", gf_field])

        if inputs.advanced_options.split_alleles:
            args.append("--split-multi-allelic")

        return_code, _, _ = Cmd["gatk"]["VariantsToTable"][args] & TEE(retcode=None)
        if return_code:
            self.error("GATK VariantsToTable failed.")

        outputs.tsv = variants_table
        outputs.species = inputs.vcf.output.species
        outputs.build = inputs.vcf.output.build

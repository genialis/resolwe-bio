"""Run variant filtration using VQSR tool."""

import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


class VariantFiltrationVqsr(Process):
    """Filter WGS variants using Variant Quality Score Recalibration (VQSR) procedure."""

    slug = "gatk-vqsr"
    name = "GATK filter variants (VQSR)"
    category = "GATK"
    process_type = "data:variants:vcf:vqsr"
    version = "1.2.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
        },
    }
    data_name = "VQSR filtered variants"

    class Input:
        """Input fields for VariantFiltrationVqsr."""

        vcf = DataField("variants:vcf", label="Input data (VCF)")

        class ResourceFiles:
            """Resource files options."""

            dbsnp = DataField("variants:vcf", label="dbSNP file")

            mills = DataField(
                "variants:vcf",
                label="Mills and 1000G gold standard indels",
                required=False,
            )

            axiom_poly = DataField(
                "variants:vcf",
                label="1000G Axiom genotype data",
                required=False,
            )

            hapmap = DataField(
                "variants:vcf",
                label="HapMap variants",
                required=False,
            )

            omni = DataField(
                "variants:vcf",
                label="1000G Omni variants",
                required=False,
            )

            thousand_genomes = DataField(
                "variants:vcf",
                label="1000G high confidence SNPs",
                required=False,
            )

        class AdvancedOptions:
            """Advanced options."""

            use_as_anno = BooleanField(
                label="--use-allele-specific-annotations", default=False
            )

            indel_anno_fields = ListField(
                StringField(),
                label="Annotation fields (INDEL filtering)",
                default=[
                    "FS",
                    "ReadPosRankSum",
                    "MQRankSum",
                    "QD",
                    "SOR",
                    "DP",
                ],
            )

            snp_anno_fields = ListField(
                StringField(),
                label="Annotation fields (SNP filtering)",
                default=[
                    "QD",
                    "MQRankSum",
                    "ReadPosRankSum",
                    "FS",
                    "MQ",
                    "SOR",
                    "DP",
                ],
            )

            indel_filter_level = FloatField(
                label="--truth-sensitivity-filter-level (INDELs)", default=99.0
            )

            snp_filter_level = FloatField(
                label="--truth-sensitivity-filter-level (SNPs)", default=99.7
            )

            max_gaussians_indels = IntegerField(
                label="--max-gaussians (INDELs)",
                default=4,
                description="This parameter determines the maximum number "
                "of Gaussians that should be used when building a positive "
                "model using the variational Bayes algorithm. This parameter "
                "sets the expected number of clusters in modeling. If a "
                "dataset gives fewer distinct clusters, e.g. as can happen "
                "for smaller data, then the tool will tell you there is "
                "insufficient data with a No data found error message. "
                "In this case, try decrementing the --max-gaussians value.",
            )

            max_gaussians_snps = IntegerField(
                label="--max-gaussians (SNPs)",
                default=6,
                description="This parameter determines the maximum number "
                "of Gaussians that should be used when building a positive "
                "model using the variational Bayes algorithm. This parameter "
                "sets the expected number of clusters in modeling. If a "
                "dataset gives fewer distinct clusters, e.g. as can happen "
                "for smaller data, then the tool will tell you there is "
                "insufficient data with a No data found error message. "
                "In this case, try decrementing the --max-gaussians value.",
            )

        resource_files = GroupField(
            ResourceFiles,
            label="Resource files",
        )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for VariantFiltrationVqsr."""

        vcf = FileField(label="GVCF file")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        variants = "snp.recalibrated.vcf"
        variants_gz = variants + ".gz"
        variants_index = variants_gz + ".tbi"

        species = inputs.vcf.output.species
        build = inputs.vcf.output.build

        excesshet_vcf = "cohort_excesshet.vcf.gz"
        sites_only_vcf = "cohort_sitesonly.vcf.gz"
        tmp_indel_recal_vcf = "indel.recalibrated.vcf"

        indels_recal = "cohort_indels.recal"
        snps_recal = "cohort_snps.recal"
        indels_tranches = "cohort_indels.tranches"
        snps_tranches = "cohort_snps.tranches"

        # Hard-filter a large cohort callset on ExcessHet
        variant_filtration = [
            "-V",
            inputs.vcf.output.vcf.path,
            "--filter-expression",
            "ExcessHet > 54.69",
            "--filter-name",
            "ExcessHet",
            "-O",
            excesshet_vcf,
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["VariantFiltration"][variant_filtration] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK VariantFiltration tool failed.")

        # Create sites-only VCF
        make_sites_only_args = [
            "-I",
            excesshet_vcf,
            "-O",
            sites_only_vcf,
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["MakeSitesOnlyVcf"][make_sites_only_args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK MakeSitesOnlyVcf tool failed.")

        indel_recalibration_tranche_values = [
            "100.0",
            "99.95",
            "99.9",
            "99.5",
            "99.0",
            "97.0",
            "96.0",
            "95.0",
            "94.0",
            "93.5",
            "93.0",
            "92.0",
            "91.0",
            "90.0",
        ]

        snp_recalibration_tranche_values = [
            "100.0",
            "99.95",
            "99.9",
            "99.8",
            "99.6",
            "99.5",
            "99.4",
            "99.3",
            "99.0",
            "98.0",
            "97.0",
            "90.0",
        ]

        args_indels = [
            "-V",
            sites_only_vcf,
            "--trust-all-polymorphic",
            "-mode",
            "INDEL",
            "--max-gaussians",
            inputs.advanced_options.max_gaussians_indels,
            "-O",
            indels_recal,
            "--tranches-file",
            indels_tranches,
            "--resource:dbsnp,known=true,training=false,truth=false,prior=2",
            inputs.resource_files.dbsnp.output.vcf.path,
            "--TMP_DIR",
            TMPDIR,
        ]

        if inputs.resource_files.mills:
            args_indels.extend(
                [
                    "--resource:mills,known=false,training=true,truth=true,prior=12",
                    inputs.resource_files.mills.output.vcf.path,
                ]
            )

        if inputs.resource_files.axiom_poly:
            args_indels.extend(
                [
                    "--resource:axiomPoly,known=false,training=true,truth=false,prior=10",
                    inputs.resource_files.axiom_poly.output.vcf.path,
                ]
            )

        for tr_val in indel_recalibration_tranche_values:
            args_indels.extend(["-tranche", tr_val])

        for anno_val in inputs.advanced_options.indel_anno_fields:
            args_indels.extend(["-an", anno_val])

        if inputs.advanced_options.use_as_anno:
            args_indels.append("--use-allele-specific-annotations")

        return_code, _, _ = Cmd["gatk"]["VariantRecalibrator"][args_indels] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK VariantRecalibrator (INDELs) tool failed.")

        args_snps = [
            "-V",
            sites_only_vcf,
            "--trust-all-polymorphic",
            "-mode",
            "SNP",
            "--max-gaussians",
            inputs.advanced_options.max_gaussians_snps,
            "-O",
            snps_recal,
            "--tranches-file",
            snps_tranches,
            "-resource:dbsnp,known=true,training=false,truth=false,prior=7",
            inputs.resource_files.dbsnp.output.vcf.path,
            "--TMP_DIR",
            TMPDIR,
        ]

        if inputs.resource_files.hapmap:
            args_snps.extend(
                [
                    "--resource:hapmap,known=false,training=true,truth=true,prior=15",
                    inputs.resource_files.hapmap.output.vcf.path,
                ]
            )

        if inputs.resource_files.omni:
            args_snps.extend(
                [
                    "-resource:omni,known=false,training=true,truth=true,prior=12",
                    inputs.resource_files.omni.output.vcf.path,
                ]
            )

        if inputs.resource_files.thousand_genomes:
            args_snps.extend(
                [
                    "-resource:1000G,known=false,training=true,truth=false,prior=10",
                    inputs.resource_files.thousand_genomes.output.vcf.path,
                ]
            )

        for tr_val in snp_recalibration_tranche_values:
            args_snps.extend(["-tranche", tr_val])

        for anno_val in inputs.advanced_options.snp_anno_fields:
            args_snps.extend(["-an", anno_val])

        if inputs.advanced_options.use_as_anno:
            args_snps.append("--use-allele-specific-annotations")

        return_code, _, _ = Cmd["gatk"]["VariantRecalibrator"][args_snps] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK VariantRecalibrator (SNPs) tool failed.")

        # ApplyVQSR
        apply_vqsr_indels = [
            "-O",
            tmp_indel_recal_vcf,
            "-V",
            excesshet_vcf,
            "--recal-file",
            indels_recal,
            "--tranches-file",
            indels_tranches,
            "--truth-sensitivity-filter-level",
            inputs.advanced_options.indel_filter_level,
            "--create-output-variant-index",
            "true",
            "-mode",
            "INDEL",
            "--TMP_DIR",
            TMPDIR,
        ]

        if inputs.advanced_options.use_as_anno:
            apply_vqsr_indels.append("--use-allele-specific-annotations")

        return_code, _, _ = Cmd["gatk"]["ApplyVQSR"][apply_vqsr_indels] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK ApplyVQSR (INDEL) tool failed.")

        apply_vqsr_snp = [
            "-O",
            variants,
            "-V",
            tmp_indel_recal_vcf,
            "--recal-file",
            snps_recal,
            "--tranches-file",
            snps_tranches,
            "--truth-sensitivity-filter-level",
            inputs.advanced_options.snp_filter_level,
            "--create-output-variant-index",
            "false",
            "-mode",
            "SNP",
            "--TMP_DIR",
            TMPDIR,
        ]

        if inputs.advanced_options.use_as_anno:
            apply_vqsr_indels.append("--use-allele-specific-annotations")

        return_code, _, _ = Cmd["gatk"]["ApplyVQSR"][apply_vqsr_snp] & TEE(retcode=None)
        if return_code:
            self.error("GATK ApplyVQSR (SNPs) tool failed.")

        # Compress and index the output variants file
        (Cmd["bgzip"]["-c", variants] > variants_gz)()
        Cmd["tabix"]["-p", "vcf", variants_gz]()

        outputs.vcf = variants_gz
        outputs.tbi = variants_index
        outputs.species = species
        outputs.build = build

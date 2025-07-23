"""Run GATK genotype refinement."""

import os

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class GatkGenotypeRefinement(Process):
    """Run GATK Genotype Refinement.

    The goal of the Genotype Refinement workflow is to use
    additional data to improve the accuracy of genotype calls
    and to filter genotype calls that are not reliable enough
    for downstream analysis. In this sense it serves as an
    optional extension of the variant calling workflow, intended
    for researchers whose work requires high-quality identification
    of individual genotypes.
    For additional information, please see
    [manual page](https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants)
    """

    slug = "gatk-refine-variants"
    name = "GATK refine variants"
    category = "GATK"
    process_type = "data:variants:vcf:refinevariants"
    version = "1.1.2"
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
    data_name = "Refined variants"
    persistence = Persistence.CACHED

    class Input:
        """Input fields for GatkGenotypeRefinement."""

        vcf = DataField(
            "variants:vcf", label="The main input, as produced in the GATK VQSR process"
        )
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

        vcf_pop = DataField(
            "variants:vcf",
            label="Population-level variant set (VCF)",
            required=False,
        )

    class Output:
        """Output fields for GatkGenotypeRefinement."""

        vcf = FileField(label="Refined multi-sample vcf")
        tbi = FileField(label="Tabix index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        variants_cgp = "variants_recal_genotype_posteriors.vcf.gz"
        variants_refined = "variants.refined.vcf.gz"
        variants_refined_index = variants_refined + ".tbi"

        if inputs.vcf_pop:
            if (
                inputs.vcf_pop.output.species
                != inputs.vcf.output.species
                != inputs.ref_seq.output.species
            ):
                self.error(
                    "Species information for the input VCF files and reference sequence don't match."
                )
            elif (
                inputs.vcf_pop.output.build
                != inputs.vcf.output.build
                != inputs.ref_seq.output.build
            ):
                self.error(
                    "Build information for the input VCF files and reference sequence don't match."
                )
        else:
            if inputs.vcf.output.species != inputs.ref_seq.output.species:
                self.error(
                    "Species information for the input VCF file and reference sequence don't match."
                )
            elif inputs.vcf.output.build != inputs.ref_seq.output.build:
                self.error(
                    "Build information for the input VCF file and reference sequence don't match."
                )

        args_cgp = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-V",
            inputs.vcf.output.vcf.path,
            "-O",
            variants_cgp,
            "--create-output-variant-index",
            "--tmp-dir",
            TMPDIR,
        ]

        if inputs.vcf_pop:
            args_cgp.extend(["--supporting", inputs.vcf_pop.output.vcf.path])

        return_code, _, _ = Cmd["gatk"]["CalculateGenotypePosteriors"][args_cgp] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK CalculateGenotypePosteriors failed.")

        args_filtration = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-V",
            variants_cgp,
            "-O",
            variants_refined,
            "--genotype-filter-name",
            "lowGQ",
            "--create-output-variant-index",
            "--genotype-filter-expression",
            "GQ < 20",
            "--tmp-dir",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["VariantFiltration"][args_filtration] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GATK VariantFiltration failed.")

        outputs.vcf = variants_refined
        outputs.tbi = variants_refined_index
        outputs.species = inputs.ref_seq.output.species
        outputs.build = inputs.ref_seq.output.build

"""Run Ensembl-VEP."""

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class EnsemblVep(Process):
    """Run Ensembl-VEP.

    VEP (Variant Effect Predictor) determines the effect of your variants
    (SNPs, insertions, deletions, CNVs or structural variants) on genes,
    transcripts, and protein sequence, as well as regulatory regions.

    This process accepts VCF file and VEP cache directory to produce
    VCF file with annotated variants, its index and summary of the procces.
    """

    slug = "ensembl-vep"
    name = "Ensembl Variant Effect Predictor"
    category = "VEP"
    process_type = "data:variants:vcf:vep"
    version = "1.1.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.2.0"}
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
            "storage": 200,
        },
    }
    data_name = '{{vcf.file|default("?") }}'

    class Input:
        """Input fields for EnsemblVep."""

        vcf = DataField(
            "variants:vcf",
            label="Input VCF file",
        )
        cache = DataField("vep:cache", label="Cache directory for Ensembl-VEP")

        n_forks = IntegerField(
            label="Number of forks",
            default=2,
            description="Using forking enables VEP to run multiple parallel threads, with each "
            "thread processing a subset of your input. Forking can dramatically improve runtime.",
        )

    class Output:
        """Output fields for EnsemblVep."""

        vcf = FileField(label="Annotated VCF file")
        tbi = FileField(label="Tabix index")
        summary = FileHtmlField(label="Summary of the analysis")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        annotated_vcf = "snp_annotated.vcf"
        annotated_vcf_gz = annotated_vcf + ".gz"
        annotated_vcf_index = annotated_vcf_gz + ".tbi"
        summary_html = annotated_vcf + "_summary.html"

        if inputs.vcf.output.species != inputs.cache.output.species:
            self.error(
                "Species information for the input VCF file and cache don't match."
            )
        elif inputs.vcf.output.build != inputs.cache.output.build:
            self.error(
                "Build information for the input VCF file and cache don't match."
            )

        if inputs.cache.output.release != "104":
            self.warning(
                "The current version of Ensembl-VEP is 104. It is recommended that cache version is also 104"
            )

        args = [
            "--cache",
            "--dir",
            inputs.cache.output.cache.path,
            "-v",
            "--format",
            "vcf",
            "--everything",
            "-i",
            inputs.vcf.output.vcf.path,
            "-o",
            annotated_vcf,
            "--fork",
            min(inputs.n_forks, self.requirements.resources.cores),
        ]

        return_code, _, _ = Cmd["vep"][args] & TEE(retcode=None)
        if return_code:
            self.error("Ensembl-VEP failed.")

        (Cmd["bgzip"]["-c", annotated_vcf] > annotated_vcf_gz)()
        Cmd["tabix"]["-p", "vcf", annotated_vcf_gz]()

        outputs.vcf = annotated_vcf_gz
        outputs.tbi = annotated_vcf_index
        outputs.summary = summary_html
        outputs.species = inputs.vcf.output.species
        outputs.build = inputs.vcf.output.build

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


def check_version(stdout):
    """Check version of Ensembl-VEP.

    Example of the first part of an output from the command `vep --help`:
    #----------------------------------#
    # ENSEMBL VARIANT EFFECT PREDICTOR #
    #----------------------------------#

    Versions:
    ensembl              : 104.1af1dce
    ensembl-funcgen      : 104.59ae779
    ensembl-io           : 104.1d3bb6e
    ensembl-variation    : 104.6154f8b
    ensembl-vep          : 104.3

    Help: dev@ensembl.org , helpdesk@ensembl.org
    Twitter: @ensembl
    """
    vep_version = int(
        float(
            next(
                (line for line in stdout.split("\n") if "ensembl-vep" in line)
            ).split()[2]
        )
    )

    return vep_version


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
    version = "2.1.0"
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
    data_name = "Annotated variants (VEP)"

    class Input:
        """Input fields for EnsemblVep."""

        vcf = DataField(
            "variants:vcf",
            label="Input VCF file",
        )
        cache = DataField("vep:cache", label="Cache directory for Ensembl-VEP")
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")

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

        _, stdout, _ = Cmd["vep"]["--help"] & TEE(retcode=None)
        vep_version = check_version(stdout=stdout)

        if inputs.cache.output.release != vep_version:
            self.warning(
                f"The current version of Ensembl-VEP is {vep_version}. "
                f"It is recommended that cache version is also {vep_version}."
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
            "--offline",
            "--fasta",
            inputs.ref_seq.output.fasta.path,
        ]

        return_code, stdout, stderr = Cmd["vep"][args] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("Ensembl-VEP failed.")

        (Cmd["bgzip"]["-c", annotated_vcf] > annotated_vcf_gz)()
        Cmd["tabix"]["-p", "vcf", annotated_vcf_gz]()

        outputs.vcf = annotated_vcf_gz
        outputs.tbi = annotated_vcf_index
        outputs.summary = summary_html
        outputs.species = inputs.vcf.output.species
        outputs.build = inputs.vcf.output.build

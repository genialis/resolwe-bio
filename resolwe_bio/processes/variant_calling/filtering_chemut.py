"""Filtering and annotation of Variant Calling (CheMut)."""
import os
from pathlib import Path

from plumbum import TEE
from pysam import VariantFile

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    IntegerField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def append_sample_info(vcf_file, summary, warning, error):
    """Extract reference (FASTA) and sample names from the VCF file."""
    try:
        vcf = VariantFile(vcf_file)
    except (OSError, ValueError) as err:
        proc_error = "Input VCF file does not exist or could not be correctly opened."
        print(err)
        error(proc_error)

    vcf_header = vcf.header
    header_records = {record.key: record.value for record in vcf_header.records}

    with open(summary, "a") as out_file:
        try:
            fasta_name = os.path.basename(header_records["reference"])
        except KeyError:
            fasta_name = ""
            warning(
                "Reference sequence (FASTA) name could not be recognized from the VCF header."
            )

        out_file.write("\nReference (genome) sequence:\n{}\n".format(fasta_name))
        out_file.write("\nSamples:\n{}".format("\n".join(list(vcf_header.samples))))


class FilteringCheMut(Process):
    """
    Filtering and annotation of Variant Calling (CheMut).

    Filtering and annotation of Variant Calling data - Chemical
    mutagenesis in _Dictyostelium discoideum_.
    """

    slug = "filtering-chemut"
    name = "Variant filtering (CheMut)"
    category = "Other"
    process_type = "data:variants:vcf:filtering"
    version = "1.7.0"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
    }
    data_name = "Filtered variants ({{ analysis_type }})"

    class Input:
        """Input fields for FilteringCheMut."""

        variants = DataField(
            data_type="variants:vcf",
            label="Variants file (VCF)",
        )
        analysis_type = StringField(
            label="Analysis type",
            description="Choice of the analysis type. Use 'SNV' or 'INDEL' options. "
            "Choose options SNV_CHR2 or INDEL_CHR2 to run the GATK analysis "
            "only on the diploid portion of CHR2 (-ploidy 2 -L chr2:2263132-3015703).",
            choices=[
                ("snv", "SNV"),
                ("indel", "INDEL"),
                ("snv_chr2", "SNV_CHR2"),
                ("indel_chr2", "INDEL_CHR2"),
            ],
            default="snv",
        )
        parental_strain = StringField(
            label="Parental strain prefix",
            default="parental",
        )
        mutant_strain = StringField(
            label="Mutant strain prefix",
            default="mut",
        )
        read_depth = IntegerField(
            label="Read Depth Cutoff",
            default=5,
        )

    class Output:
        """Output fields for FilteringCheMut."""

        summary = FileField(
            label="Summary", description="Summarize the input parameters and results."
        )
        vcf = FileField(
            label="Variants",
            description="A genome VCF file of variants that passed the filters.",
        )
        tbi = FileField(label="Tabix index")
        variants_filtered = FileField(
            label="Variants filtered",
            required=False,
            description="A data frame of variants that passed the filters.",
        )
        variants_filtered_alt = FileField(
            label="Variants filtered (multiple alt. alleles)",
            required=False,
            description="A data frame of variants that contain more than two alternative "
            "alleles. These variants are likely to be false positives.",
        )
        gene_list_all = FileField(
            label="Gene list (all)",
            required=False,
            description="Genes that are mutated at least once.",
        )
        gene_list_top = FileField(
            label="Gene list (top)",
            required=False,
            description="Genes that are mutated at least twice.",
        )
        mut_chr = FileField(
            label="Mutations (by chr)",
            required=False,
            description="List mutations in individual chromosomes.",
        )
        mut_strain = FileField(
            label="Mutations (by strain)",
            required=False,
            description="List mutations in individual strains.",
        )
        strain_by_gene = FileField(
            label="Strain (by gene)",
            required=False,
            description="List mutants that carry mutations in individual genes.",
        )
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        base_path = Path(inputs.variants.output.vcf.path)
        assert base_path.name.endswith(".vcf.gz")
        vcf_file = base_path.stem

        (Cmd["bgzip"]["-dc"][inputs.variants.output.vcf.path] > vcf_file)()

        r_input = (
            "library(chemut); "
            f"{inputs.analysis_type}("
            f"input_file = '{vcf_file}', "
            f"parental_strain = '{inputs.parental_strain}', "
            f"mutant_strain = '{inputs.mutant_strain}', "
            f"read_depth = {inputs.read_depth})"
        )

        return_code, _, stderr = Cmd["Rscript"]["-e"][r_input] & TEE(retcode=None)
        if return_code:
            print(stderr)
            self.error(f"Error while running the script {inputs.analysis_type}.R")

        output_dir = Path(f"{vcf_file}_{inputs.read_depth}")
        if not output_dir.exists():
            output_dir.mkdir()

        append_sample_info(
            vcf_file=vcf_file,
            summary=str(output_dir / "summary.txt"),
            warning=self.warning,
            error=self.error,
        )

        outputs.summary = str(output_dir / "summary.txt")
        outputs.species = inputs.variants.output.species
        outputs.build = inputs.variants.output.build

        if (output_dir / "variants.vcf").exists():

            variants_gz = str(output_dir / "variants.vcf.gz")
            (Cmd["bgzip"]["-c", str(output_dir / "variants.vcf")] > variants_gz)()

            Cmd["tabix"]["-p", "vcf", variants_gz]()
            outputs.vcf = variants_gz
            outputs.tbi = variants_gz + ".tbi"
        else:
            self.error("No variants have passed the filters. VCF file was not created.")

        if (output_dir / "variants_filtered.txt").exists():
            outputs.variants_filtered = str(output_dir / "variants_filtered.txt")

        if (output_dir / "variant_mult_alt.txt").exists():
            outputs.variants_filtered_alt = str(output_dir / "variant_mult_alt.txt")

        if (output_dir / "gene_list_all.txt").exists():
            outputs.gene_list_all = str(output_dir / "gene_list_all.txt")

        if (output_dir / "gene_list_top.txt").exists():
            outputs.gene_list_top = str(output_dir / "gene_list_top.txt")

        if (output_dir / "mutations_by_chr.txt").exists():
            outputs.mut_chr = str(output_dir / "mutations_by_chr.txt")

        if (output_dir / "mutations_by_strain.txt").exists():
            outputs.mut_strain = str(output_dir / "mutations_by_strain.txt")

        if (output_dir / "strain_by_gene.txt").exists():
            outputs.strain_by_gene = str(output_dir / "strain_by_gene.txt")

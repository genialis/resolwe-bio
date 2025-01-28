"""Variant annotation process."""

import gzip
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    Persistence,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio

VCF_FIELDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "ANN",
    "CLNDN",
    "CLNSIG",
]

ANN_FIELDS = [
    "Allele",
    "Annotation",
    "Annotation_Impact",
    "Gene_Name",
    "Feature_ID",
    "HGVS.p",
]


ANN_COLUMNS = [
    "Allele",
    "Annotation",
    "Annotation_Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS.c",
    "HGVS.p",
    "cDNA.pos/cDNA.length",
    "CDS.pos/CDS.length",
    "AA.pos/AA.length",
    "Distance",
    "ERRORS/WARNINGS/INFO",
]


def extract_ids(id_column):
    """Extract SNP and ClinVar IDs from ID column."""
    id_column = id_column.str.replace(".", "", regex=False)
    ids_expanded = id_column.str.split(";", expand=True)

    def extract_snp_ids(x):
        return ",".join(x[x.str.startswith("rs", na=False)].dropna())

    def extract_clinvar_ids(x):
        return ",".join(x[~x.str.startswith("rs", na=False)].dropna())

    snpids = ids_expanded.apply(extract_snp_ids, axis=1)
    clinvarids = ids_expanded.apply(extract_clinvar_ids, axis=1)

    return pd.DataFrame({"SNPID": snpids, "CLINVARID": clinvarids})


def prepare_canonical(geneset):
    """Prepare gene set for further analysis."""
    transcripts = []
    with gzip.open(geneset, "rb") as file_in:
        transcripts = [transcript.decode().rstrip() for transcript in file_in]
    return transcripts


def compare_ref_alt(ref, alt):
    """Comparison of reference and alternative alleles."""
    symbolic_mask = alt.str.startswith("<") & alt.str.endswith(">")
    snp_mask = (ref.str.len() == 1) & (alt.str.len() == 1) & ~symbolic_mask
    mnp_mask = (ref.str.len() == alt.str.len()) & (ref.str.len() > 1) & ~symbolic_mask

    result = pd.Series(["INDEL"] * len(ref))

    result[symbolic_mask] = "SYMBOLIC"
    result[snp_mask] = "SNP"
    result[mnp_mask] = "MNP"

    return result


def infer_variant_type(ref, alt):
    """Inference of variant type.

    This code flattens the variants into a 1D array and assigns
    index for each VCF entry (variant_indices).
    Multiallelic variants will therefore have duplicated indices.
    Then, it groups the entries by the index and determines the variant type for each group.
    """

    def determine_mixed_variants(vcf_entry):
        """Determine if variant call is of a mixed type.

        If there are multiple different variant types present in a multiallelic call,
        the variant is considered mixed.
        """
        unique_types = set(vcf_entry)
        return "MIXED" if len(unique_types) > 1 else unique_types.pop()

    alt_split = alt.str.split(",")
    refs = ref.repeat(alt_split.str.len())
    alts = np.concatenate(alt_split.values)
    variant_indices = np.repeat(alt_split.index, alt_split.str.len())

    types = compare_ref_alt(ref=pd.Series(refs.values), alt=pd.Series(alts))
    inferred_types = types.groupby(variant_indices).agg(determine_mixed_variants)

    return inferred_types


def process_variant(variant_df, canonical_transcripts, genome_assembly):
    """Process a single VCF entry."""

    variant_dict = {
        "species": "Homo sapiens",
        "genome_assembly": genome_assembly,
        "chromosome": variant_df["CHROM"].iloc[0],
        "position": int(variant_df["POS"].iloc[0]),
        "reference": variant_df["REF"].iloc[0],
        "alternative": variant_df["ALT"].iloc[0],
        "type": variant_df["TYPE"].iloc[0],
        "clinical_diagnosis": variant_df["CLNDN"].iloc[0],
        "clinical_significance": variant_df["CLNSIG"].iloc[0],
        "dbsnp_id": variant_df["SNPID"].iloc[0],
        "clinvar_id": variant_df["CLINVARID"].iloc[0],
        "transcripts": [
            {
                "annotation": row.Annotation,
                "annotation_impact": row.Annotation_Impact,
                "gene": row.Gene_Name,
                "protein_impact": row.HGVS_p,
                "transcript_id": row.Feature_ID.split(".")[0],
                "canonical": row.Feature_ID.split(".")[0] in canonical_transcripts,
            }
            for row in variant_df.itertuples(index=False, name="Pandas")
        ],
    }

    return variant_dict


def create_annotation_list(variants_table, canonical_transcripts, genome_assembly):
    """Create annotation list."""
    variants_table = variants_table.astype(str)
    column_names = {name: name.replace(".", "_") for name in variants_table.columns}
    variants_table.columns = [column_names[col] for col in variants_table.columns]
    variants_table = variants_table.replace("nan", "")

    variant_group = variants_table.groupby(
        ["CHROM", "POS", "REF", "ALT", "TYPE", "CLNDN", "CLNSIG", "SNPID", "CLINVARID"]
    )

    annotations_list = Parallel(n_jobs=-1)(
        delayed(process_variant)(
            variant_df=variant_df,
            canonical_transcripts=canonical_transcripts,
            genome_assembly=genome_assembly,
        )
        for _, variant_df in variant_group
    )

    return annotations_list


def prepare_variants_table(variants_table, ann_fields, vcf_fields):
    """Prepare variants table."""
    variants_table = pd.read_csv(
        variants_table,
        sep="\t",
        dtype={
            "CHROM": str,
            "POS": int,
            "REF": str,
            "ALT": str,
            "QUAL": float,
            "ANN": str,
            "CLNDN": str,
            "CLNSIG": str,
        },
    )
    # First split each line to multiple lines, since every variant has
    # info for multiple transcripts.
    variants_table = variants_table.drop("ANN", axis=1).join(
        variants_table["ANN"]
        .str.split(",", expand=True)
        .stack()
        .reset_index(level=1, drop=True)
        .rename("ANN")
    )
    # Split ANN column to multiple columns
    variants_table[ANN_COLUMNS] = variants_table["ANN"].str.split("|", expand=True)
    # Only use subset of the original dataframe
    variants_table = variants_table[vcf_fields + ann_fields]
    variants_table = variants_table.drop(["ANN", "QUAL"], axis=1)
    variants_table.reset_index(inplace=True, drop=True)

    variants_table[["SNPID", "CLINVARID"]] = extract_ids(id_column=variants_table["ID"])
    variants_table["TYPE"] = infer_variant_type(
        ref=variants_table["REF"], alt=variants_table["ALT"]
    )

    variants_table = variants_table.drop_duplicates()

    return variants_table


class VariantAnnotation(ProcessBio):
    """Variant annotation process."""

    slug = "variant-annotation"
    name = "Variant annotation"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/snpeff:3.2.0",
            },
        },
        "resources": {
            "cores": 8,
            "memory": 8196,
        },
    }
    data_name = "Variant annotation"
    version = "1.2.0"
    process_type = "data:annotation"
    category = "Annotation"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.TEMP

    description = "Variant annotation using SnpEff and ClinVar."

    class Input:
        """Input fields."""

        dbsnp = DataField(
            data_type="variants:vcf",
            label="dbSNP VCF file",
            description="Database of known polymorphic sites.",
            required=True,
        )
        clinvar = DataField(
            data_type="variants:vcf",
            label="ClinVar VCF file",
            description="[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a "
            "freely available, public archive of human genetic variants and "
            "interpretations of their significance to disease.",
            required=True,
        )
        canonical_transcripts = DataField(
            data_type="geneset",
            label="Canonical transcripts",
            description="List of canonical transcripts for annotation.",
            required=True,
        )
        database = StringField(
            label="snpEff database",
            default="GRCh38.109",
            choices=[
                ("GRCh37.75", "GRCh37.75"),
                ("GRCh38.99", "GRCh38.99"),
                ("GRCh38.109", "GRCh38.109"),
            ],
        )
        select_all = BooleanField(
            label="Select all variants for annotation",
            description="Select all variants for annotation, or just those that do not exist in the database.",
            default=False,
        )
        genome_assembly = StringField(
            label="Genome assembly",
            description="Genome assembly",
            choices=[("GRCh37", "GRCh37"), ("GRCh38", "GRCh38")],
            default="GRCh38",
        )

        class Advanced:
            """Advanced options."""

            save_vcf = BooleanField(
                label="Save annotated VCF",
                default=False,
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields."""

        vcf = FileField(label="Annotated VCF file", required=False)
        tbi = FileField(label="Tabix index", required=False)
        genes = FileField(label="SnpEff genes")
        summary = FileField(label="SnpEff summary")

    def run(self, inputs, outputs):
        """Run the process."""

        snpeff_variants = "snpeff_variants.vcf"
        dbsnp_variants = "dbsnp_variants.vcf"
        annotated_variants = "annotated_variants.vcf"

        filters = {
            "genome_assembly": inputs.genome_assembly,
            "__fields": ["chromosome", "position", "reference", "alternative"],
        }
        if not inputs.select_all:
            filters["annotation__isnull"] = True

        variants = [variant for variant in self.variant.iterate(**filters)]

        if not variants:
            self.error("No variants found in the database for the selected inputs.")

        canonical_transcripts = prepare_canonical(
            geneset=inputs.canonical_transcripts.output.geneset.path
        )

        vcf_data = [
            {
                "CHROM": variant.chromosome,
                "POS": variant.position,
                "ID": ".",
                "REF": variant.reference,
                "ALT": variant.alternative,
                "QUAL": ".",
                "FILTER": ".",
                "INFO": ".",
                "FORMAT": ".",
                "SAMPLENAME1": ".",
            }
            for variant in variants
        ]

        vcf_df = pd.DataFrame(vcf_data)
        vcf_df = vcf_df.sort_values(["CHROM", "POS"], ascending=[True, True])
        vcf_fn = "variants.vcf"
        vcf_df.to_csv("variants.vcf", sep="\t", index=False, header=False)
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLENAME1"
        with open(vcf_fn, "r") as original:
            data = original.read()
        with open(vcf_fn, "w") as modified:
            modified.write(header + "\n" + data)

        self.progress(0.1)

        args_snpeff = [
            inputs.database,
            vcf_fn,
            "-no-downstream",
            "-no-upstream",
            "-no-intergenic",
        ]
        (Cmd["snpEff"][args_snpeff] > snpeff_variants)()

        self.progress(0.3)

        if not inputs.dbsnp.output.build.startswith(inputs.database[:6]):
            self.error(
                "Genome build for the DBSNP file and used database "
                "should be the same. DBSNP file is based on "
                f"{inputs.dbsnp.output.build}, while snpEff database "
                f"is based on {inputs.database[:6]}."
            )
        args_annotation = [
            "annotate",
            inputs.dbsnp.output.vcf.path,
            snpeff_variants,
        ]

        (Cmd["SnpSift"][args_annotation] > dbsnp_variants)()

        if not inputs.clinvar.output.build.startswith(inputs.database[:6]):
            self.error(
                "Genome build for the ClinVar file and used database "
                "should be the same. ClinVar file is based on "
                f"{inputs.clinvar.output.build}, while snpEff database "
                f"is based on {inputs.database[:6]}."
            )
        args_annotation = [
            "annotate",
            inputs.clinvar.output.vcf.path,
            dbsnp_variants,
        ]

        (Cmd["SnpSift"][args_annotation] > annotated_variants)()

        with open(annotated_variants, "r+") as file:
            content = file.read()
            file.seek(0)
            file.write("##fileformat=VCFv4.2\n" + content)

        gzipped_variants = f"{annotated_variants}.gz"
        (Cmd["bgzip"]["-c", annotated_variants] > gzipped_variants)()
        (Cmd["tabix"]["-p", "vcf", gzipped_variants])()

        if inputs.advanced.save_vcf:
            outputs.vcf = gzipped_variants
            outputs.tbi = f"{gzipped_variants}.tbi"

        outputs.genes = "snpEff_genes.txt"
        outputs.summary = "snpEff_summary.html"

        self.progress(0.5)

        TMPDIR = os.environ.get("TMPDIR")
        variants_table = "variants_table.tsv"

        args = [
            "-V",
            gzipped_variants,
            "-O",
            variants_table,
            "--tmp-dir",
            TMPDIR,
        ]

        for field in VCF_FIELDS:
            args.extend(["-F", field])

        args.append("--split-multi-allelic")
        args.append("--show-filtered")

        return_code, _, _ = Cmd["gatk"]["VariantsToTable"][args] & TEE(retcode=None)
        if return_code:
            self.error("GATK VariantsToTable failed.")

        self.progress(0.7)

        variants_table = prepare_variants_table(
            variants_table=variants_table, ann_fields=ANN_FIELDS, vcf_fields=VCF_FIELDS
        )

        annotations_list = create_annotation_list(
            variants_table=variants_table,
            canonical_transcripts=canonical_transcripts,
            genome_assembly=inputs.genome_assembly,
        )

        self.add_variants_annotations(annotations_list)

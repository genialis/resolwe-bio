"""Report mutations from RNA-seq Variant Calling Pipeline."""

import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    ListField,
    Persistence,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio

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

AMINOACIDS = [
    "Arg",
    "His",
    "Lys",
    "Asp",
    "Glu",
    "Ser",
    "Thr",
    "Asn",
    "Gln",
    "Gly",
    "Pro",
    "Cys",
    "Ala",
    "Val",
    "Ile",
    "Leu",
    "Met",
    "Phe",
    "Tyr",
    "Trp",
]


def get_output_table(mutations, variants_table, output_table, warning):
    """Prepare output table."""

    genes = []

    for gene in mutations:
        genes.append(gene)
        if len(mutations[gene]) > 0:
            for mutation in mutations[gene]:
                df = variants_table.loc[
                    (variants_table["HGVS.p"].str.contains(mutation))
                    & (variants_table["Gene_Name"] == gene)
                ]
                output_table = pd.concat([output_table, df], ignore_index=True, axis=0)
        else:
            df = variants_table.loc[(variants_table["Gene_Name"] == gene)]
            output_table = pd.concat([output_table, df], ignore_index=True, axis=0)

    if output_table.empty and not variants_table.empty:
        warning("No variants present for the input set of mutations.")

    return output_table, genes


def get_mutations(input_mutations, error):
    """Input mutations to dictionary."""

    mutations = defaultdict(list)
    for mutation in input_mutations:
        mutation = mutation.replace(" ", "")
        gene = mutation.split(sep=":")[0]
        if len(mutation.split(sep=":")) == 2:
            aminoacids = mutation.split(sep=":")[1].split(sep=",")
            for aminoacid in aminoacids:
                if aminoacid[:3] not in AMINOACIDS:
                    error(
                        f"The input amino acid {aminoacid[:3]} is in the wrong format "
                        "or is not among the 20 standard amino acids."
                    )
                mutations[gene].append(aminoacid)
        elif len(mutation.split(sep=":")) == 1:
            mutations[gene] = []
        else:
            error("Wrong input format for mutations.")

    return mutations


def prepare_geneset(geneset):
    """Prepare gene set for further analysis."""

    mutations = []
    with gzip.open(geneset, "rb") as file_in:
        for gene in file_in:
            mutations.append(gene.decode().rstrip())

    return mutations


def prepare_variants_table(variants_table, vcf_fields, ann_fields, gt_fields, warning):
    """Prepare variants table."""

    variants_table = pd.read_csv(
        variants_table,
        sep="\t",
        header=0,
        float_precision="round_trip",
    )
    if variants_table.empty:
        warning("There are no variants in the input VCF file.")

        for ann in ann_fields:
            variants_table[ann] = None

        variants_table.drop("ANN", axis=1, inplace=True)
    else:
        variants_table = ann_field_to_df(
            variants_table=variants_table,
            ann_fields=ann_fields,
            vcf_fields=vcf_fields,
            gt_fields=[f"SAMPLENAME1.{field}" for field in gt_fields],
        )
        # Collapse multiple transcripts of one variant to one line
        variants_table = (
            variants_table.groupby(
                vcf_fields
                + [field for field in ann_fields if field != "Feature_ID"]
                + [
                    col
                    for col in variants_table.columns
                    if col.startswith("SAMPLENAME1")
                ],
                dropna=False,
            )["Feature_ID"]
            .apply(",".join)
            .reset_index()
        )
    # Migrate sample-level filter value to the general FILTER field
    if "FT" in gt_fields:
        variants_table["FILTER"] = np.where(
            (
                ~variants_table["SAMPLENAME1.FT"].isna()
                & ~variants_table["FILTER"].str.contains("PASS")
            ),
            variants_table["FILTER"] + ";" + variants_table["SAMPLENAME1.FT"],
            variants_table["FILTER"],
        )
        variants_table["FILTER"] = np.where(
            (
                ~variants_table["SAMPLENAME1.FT"].isna()
                & variants_table["FILTER"].str.contains("PASS")
            ),
            variants_table["SAMPLENAME1.FT"],
            variants_table["FILTER"],
        )
        variants_table.drop(["SAMPLENAME1.FT"], axis=1, inplace=True)

    if "DP" in gt_fields:
        if "DP" in vcf_fields:
            variants_table.drop(["DP"], axis=1, inplace=True)

        variants_table.rename(columns={"SAMPLENAME1.DP": "DP"}, inplace=True)

    return variants_table


def ann_field_to_df(variants_table, ann_fields, vcf_fields, gt_fields=None):
    """Transform SnpEff ANN field to multiple rows and columns."""

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
    if gt_fields:
        variants_table = variants_table[vcf_fields + ann_fields + gt_fields]
    else:
        variants_table = variants_table[vcf_fields + ann_fields]
    variants_table.reset_index(inplace=True, drop=True)

    return variants_table


def get_depth(variants_table, bam):
    """Calculate depth for every positions in the variants table."""

    bases = ["Base_A", "Base_C", "Base_G", "Base_T"]

    if not variants_table.empty:
        for i in range(len(bases)):
            variants_table[bases[i]] = variants_table.apply(
                lambda row: bam.count_coverage(
                    contig=str(row["CHROM"]), start=row["POS"] - 1, stop=row["POS"]
                )[i][0],
                axis=1,
            )
        variants_table["Total_depth"] = variants_table[bases].sum(axis=1)
        variants_table["POS"] = variants_table["POS"].astype(int)
        variants_table.sort_values(by=["POS"], inplace=True)

    else:
        for col in bases + ["Total_depth"]:
            variants_table[col] = None

    return variants_table.to_csv("mutations.tsv", sep="\t")


def encode_variant_genotype(row):
    """Encode variant to genotype notation.

    Variants are given as <allele1>/<allele2>, e.g. T/T or C/T

    Encode these variants as:
        - 0/1 for heterozygous variants
        - 1/1 for homozygous variants
    """
    try:
        allele_line = row.get("genotype", np.nan)
        allele_re = r"^([ATGC*]+)/([ATGC*]+)$"
        allele1, allele2 = re.match(allele_re, allele_line).group(1, 2)
    except AttributeError:
        # AttributeError is raised when there is no match, e.g.
        # there is a string value for column "genotype" but
        # the above regex can't parse it
        print(f'Cannot encode variant from value "{allele_line}".')
        return np.nan

    if allele1 == allele2 == row["alternative"]:
        return "1/1"
    else:
        return "0/1"


class MutationsTable(ProcessBio):
    """Report mutations in a table from RNA-seq Variant Calling Pipeline.

    This process reports only mutations selected in the process input.
    For example, if you want to know, if the mutation Gly12X in the gene
    KRAS is present in the sample, use VCF file annotated with SnpEff as
    input together with the desired gene and mutation. One of the inputs
    should also be annotated dbSNP VCF file. The process also calculates
    sequencing coverage at positions of desired mutations, so one of the
    inputs should also be a BAM file from BQSR or SplitNCigarReads process.
    """

    slug = "mutations-table"
    name = "Mutations table"
    process_type = "data:mutationstable"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"},
        },
        "resources": {
            "cores": 1,
            "memory": 8196,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "WGS"
    data_name = "{{ variants|name|default('?') }}"
    version = "2.3.0"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    class Input:
        """Input fields to ReportVariants."""

        variants = DataField(
            data_type="variants:vcf:snpeff",
            label="Annotated variants",
            description="Variants annotated with SnpEff. VCF file used for "
            "annotation should only be filtered but should include the filtered "
            "variants.",
        )
        mutations = ListField(
            StringField(),
            label="Gene and its mutations",
            description="Insert the gene you are interested in, together "
            "with mutations. First enter the name of the gene and then "
            "the mutations. Seperate gene from mutations with ':' and mutations "
            "with ','. Example of an input: 'KRAS: Gly12, Gly61'. Press enter "
            "after each input (gene + mutations). NOTE: Field only accepts "
            "three character amino acid symbols.",
            disabled="geneset",
            required=False,
        )
        geneset = DataField(
            data_type="geneset",
            label="Gene set",
            description="Select a gene set with genes you are interested in. "
            "Only variants of genes in the selected gene set will be in the "
            "output.",
            disabled="mutations",
            required=False,
        )
        vcf_fields = ListField(
            StringField(),
            label="Select VCF fields",
            description="The name of a standard VCF field or an "
            "INFO field to include in the output table. "
            "The field can be any standard VCF column (e.g. CHROM, ID, QUAL) "
            "or any annotation name in the INFO field (e.g. AC, AF). "
            "Required fields are CHROM, POS, ID, REF, ALT and ANN. If your variants "
            "file was annotated with clinvar information then fields CLNDN, "
            "CLNSIG and CLNSIGCONF might be of your interest.",
            default=[
                "CHROM",
                "POS",
                "ID",
                "QUAL",
                "REF",
                "ALT",
                "DP",
                "QD",
                "FILTER",
                "ANN",
            ],
        )
        ann_fields = ListField(
            StringField(),
            label="ANN fields to use",
            description="Only use specific fields from the SnpEff ANN "
            "field. All available fields: Allele | Annotation | Annotation_Impact "
            "| Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType "
            "| Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length "
            "| AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ."
            "Fields are seperated by '|'. For more information, follow this [link]"
            "(https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files).",
            default=[
                "Allele",
                "Annotation",
                "Annotation_Impact",
                "Gene_Name",
                "Feature_ID",
                "HGVS.p",
            ],
        )
        bam = DataField(
            data_type="alignment:bam",
            label="Bam file used for coverage calculation",
            description="Output BAM file from BQSR or SplitNCigarReads should be used. BAM "
            "file should be from the same sample as the input variants file.",
        )

        class Advanced:
            """Advanced options."""

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
            show_filtered = BooleanField(
                label="Include filtered records in the output",
                default=True,
                description="Include filtered records in the output of the GATK "
                "VariantsToTable.",
            )
            gf_fields = ListField(
                StringField(),
                label="Include FORMAT/sample-level fields. Note: If you specify DP "
                "from genotype field, it will overwrite the original DP field.",
                default=[
                    "GT",
                    "GQ",
                    "AD",
                ],
            )
            data_source = StringField(
                label="Data source",
                description="Data source of the reported variants.",
                default="RNA-Seq Variant Calling Pipeline",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to ReportVariants."""

        tsv = FileField(
            label="Mutations table",
        )
        genes = ListField(StringField(), label="Input genes")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        if not inputs.mutations and not inputs.geneset:
            self.error(
                "Mutations or geneset were not specified. You must either enter desired "
                "mutations or select your geneset of interest."
            )
        if not all(
            field in inputs.vcf_fields
            for field in ["CHROM", "POS", "ID", "REF", "ALT", "ANN"]
        ):
            self.error(
                "Input VCF fields do not contain all required values. "
                "Required fields are CHROM, POS, ID, REF, ALT and ANN."
            )

        if inputs.variants.entity.id != inputs.bam.entity.id:
            self.error(
                "Sample ids of input annotated variants and input bam file do not match. "
                f"Annotated variants have sample id {inputs.variants.entity.id}, while bam "
                f"file has sample id {inputs.bam.entity.id}."
            )

        TMPDIR = os.environ.get("TMPDIR")
        variants_table = "variants_table.tsv"

        args = [
            "-V",
            inputs.variants.output.vcf.path,
            "-O",
            variants_table,
            "--tmp-dir",
            TMPDIR,
        ]

        for field in inputs.vcf_fields:
            args.extend(["-F", field])
        for gf_field in inputs.advanced.gf_fields:
            args.extend(["-GF", gf_field])
        if inputs.advanced.split_alleles:
            args.append("--split-multi-allelic")
        if inputs.advanced.show_filtered:
            args.append("--show-filtered")

        return_code, stdout, stderr = Cmd["gatk"]["VariantsToTable"][args] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("GATK VariantsToTable failed.")

        vcf_fields = inputs.vcf_fields
        vcf_fields.remove("ANN")

        variants_table = prepare_variants_table(
            variants_table=variants_table,
            vcf_fields=vcf_fields,
            ann_fields=inputs.ann_fields,
            gt_fields=inputs.advanced.gf_fields,
            warning=self.warning,
        )

        if inputs.mutations:
            mutations = get_mutations(
                input_mutations=inputs.mutations, error=self.error
            )
        elif inputs.geneset:
            geneset = prepare_geneset(inputs.geneset.output.geneset.path)

            feature_filters = {
                "source": inputs.geneset.output.source,
                "species": inputs.geneset.output.species,
                "feature_id__in": geneset,
            }

            geneset = [f.name for f in self.feature.filter(**feature_filters)]
            if len(geneset) == 0:
                self.error(
                    "Geneset is either empty or no gene IDs were mapped to gene symbols."
                )

            mutations = defaultdict(list)
            for gene in geneset:
                mutations[gene] = []

        bam = pysam.AlignmentFile(inputs.bam.output.bam.path, "rb")
        # Set up the output dataframe
        output_table = pd.DataFrame()

        output_table, genes = get_output_table(
            mutations=mutations,
            variants_table=variants_table,
            output_table=output_table,
            warning=self.warning,
        )

        get_depth(variants_table=output_table, bam=bam)
        output_table.reset_index(inplace=True, drop=True)

        output_table.to_csv("mutations.tsv", sep="\t", index=False)

        outputs.tsv = "mutations.tsv"
        outputs.genes = genes
        outputs.species = inputs.variants.output.species
        outputs.build = inputs.variants.output.build

        # Store the reported variants in the Variant Database
        output_table["species"] = inputs.variants.output.species
        output_table["genome_assembly"] = (
            "GRCh38"
            if "GRCh38" in inputs.variants.output.build
            else inputs.variants.output.build
        )

        columns_map = {
            "species": "species",
            "genome_assembly": "genome_assembly",
            "CHROM": "chromosome",
            "POS": "position",
            "REF": "reference",
            "ALT": "alternative",
        }

        optional_mapping = {
            "QUAL": "quality",
            "DP": "depth",
            "FILTER": "filter",
            "QD": "depth_norm_quality",
            "SAMPLENAME1.AD": "unfiltered_allele_depth",
            "SAMPLENAME1.GT": "genotype",
            "SAMPLENAME1.GQ": "genotype_quality",
        }

        gf_fields = [f"SAMPLENAME1.{field}" for field in inputs.advanced.gf_fields]
        requested_fields = inputs.vcf_fields + gf_fields

        if "DP" in inputs.advanced.gf_fields:
            requested_fields.append("DP")

        for key, value in optional_mapping.items():
            if key in requested_fields:
                columns_map[key] = value

        output_table = output_table.rename(columns=columns_map)

        if "genotype" in output_table.columns:
            output_table["genotype"] = output_table.apply(
                encode_variant_genotype, axis=1, result_type="reduce"
            )

        if "unfiltered_allele_depth" in output_table.columns:
            output_table["unfiltered_allele_depth"] = (
                output_table["unfiltered_allele_depth"]
                .str.split(",")
                .str[1]
                .astype(int)
            )

        # before reporting sample variants, drop columns containing variant
        # annotation data, e.g "HGVS.p". This allows duplicate variant rows
        # to be dropped.
        retain_fields = list(columns_map.values())
        output_table = output_table[retain_fields]
        output_table = output_table.drop_duplicates()
        output_table = output_table[columns_map.values()].to_dict(orient="records")
        self.add_variants(inputs.advanced.data_source, output_table)

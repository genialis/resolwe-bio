"""Report mutations from RNA-seq Variant Calling Pipeline."""
import gzip
import os
from collections import defaultdict

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


def get_output_table(mutations, reference, variants_table, output_table, genes, error):
    """Prepare output table."""

    for gene in mutations:
        genes.append(gene)
        if len(mutations[gene]) > 0:
            for mutation in mutations[gene]:
                # Mutations for specific codon can theoretically be at three
                # different sites
                positions = list(
                    set(
                        reference.loc[
                            (reference["HGVS.p"].str.contains(mutation))
                            & (reference["Gene_Name"] == gene),
                            "POS",
                        ].tolist()
                    )
                )
                if len(positions) == 0:
                    error(
                        f"There are no known variants for gene {gene} at "
                        f"amino acid {mutation}."
                    )
                if (
                    variants_table.empty
                    or (
                        variants_table.loc[
                            (variants_table["HGVS.p"].str.contains(mutation))
                            & (variants_table["Gene_Name"] == gene)
                        ]
                    ).empty
                ):
                    for position in positions:
                        # Using [0] here because because if reference table was processed with
                        # the option 'split_alleles', there will be multiple rows for one position
                        chrom = reference.loc[
                            reference["POS"] == position, "CHROM"
                        ].tolist()[0]
                        id = reference.loc[reference["POS"] == position, "ID"].tolist()[
                            0
                        ]
                        ref = reference.loc[
                            reference["POS"] == position, "REF"
                        ].tolist()[0]
                        df = pd.DataFrame(
                            {
                                "CHROM": [chrom],
                                "POS": [position],
                                "ID": [id],
                                "REF": [ref],
                                "Gene_Name": [gene],
                                "HGVS.p": [
                                    f"no mutation at position {mutation} detected"
                                ],
                            },
                        )
                        output_table = pd.concat(
                            [output_table, df], ignore_index=True, axis=0
                        )
                else:
                    df = variants_table.loc[
                        (variants_table["HGVS.p"].str.contains(mutation))
                        & (variants_table["Gene_Name"] == gene)
                    ]
                    output_table = pd.concat(
                        [output_table, df], ignore_index=True, axis=0
                    )
        else:
            positions = list(
                set(
                    reference.loc[
                        (reference["Gene_Name"] == gene),
                        "POS",
                    ].tolist()
                )
            )
            if len(positions) == 0:
                error(f"There are no known variants for gene {gene}.")
            if (
                variants_table.empty
                or (variants_table.loc[(variants_table["Gene_Name"] == gene)]).empty
            ):
                for position in positions:
                    chrom = reference.loc[
                        reference["POS"] == position, "CHROM"
                    ].tolist()[0]
                    id = reference.loc[reference["POS"] == position, "ID"].tolist()[0]
                    ref = reference.loc[reference["POS"] == position, "REF"].tolist()[0]
                    df = pd.DataFrame(
                        {
                            "CHROM": [chrom],
                            "POS": [position],
                            "ID": [id],
                            "REF": [ref],
                            "Gene_Name": [gene],
                        },
                    )
                    output_table = pd.concat(
                        [output_table, df], ignore_index=True, axis=0
                    )
            else:
                df = variants_table.loc[(variants_table["Gene_Name"] == gene)]
                output_table = pd.concat([output_table, df], ignore_index=True, axis=0)

    return output_table, genes


def check_reference(reference, error):
    """Check if reference contains all necessary VCF fields."""

    if (
        not pd.Series(["CHROM", "POS", "ID", "REF", "ANN"])
        .isin(reference.columns)
        .all()
    ):
        error(
            "Reference variants table does not contain all necessary fields. "
            "Necessary fields are CHROM, POS, ID, REF and ANN."
        )


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
    else:
        vcf_fields.remove("ANN")
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
                + [col for col in variants_table.columns if col.endswith("GT")],
                dropna=False,
            )["Feature_ID"]
            .apply(",".join)
            .reset_index()
        )

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

    for i in range(len(bases)):
        variants_table[bases[i]] = variants_table.apply(
            lambda row: bam.count_coverage(
                contig=row["CHROM"], start=row["POS"] - 1, stop=row["POS"]
            )[i][0],
            axis=1,
        )
    variants_table["Total_depth"] = variants_table[bases].sum(axis=1)
    variants_table["POS"] = variants_table["POS"].astype(int)
    variants_table.sort_values(by=["POS"], inplace=True)

    return variants_table.to_csv("mutations.tsv", sep="\t")


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
    category = "Other"
    data_name = "{{ variants|name|default('?') }}"
    version = "1.3.0"
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
            "Required fields are CHROM, POS, ID, REF and ANN. If your variants "
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
        reference = DataField(
            data_type="variantstable",
            label="Table of known annotated variants",
            description="To get the appropriate input format, first reference file (dbSNP) "
            "has to be annotated with SnpEff. Then annotated VCF file has to be processed "
            "by the process GATK VariantsToTable. When triggering the process "
            "VariantsToTable keep in mind that the process Mutations table needs fields "
            "CHROM, POS, ID, REF and ANN in the input reference file. Output table of "
            "the process GATK VariantsToTable can be used as an input for this field.",
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
                label="Include FORMAT/sample-level fields",
                default=[
                    "GT",
                ],
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
            field in inputs.vcf_fields for field in ["CHROM", "POS", "ID", "REF", "ANN"]
        ):
            self.error(
                "Input VCF fields do not contain all required values. "
                "Required fields are CHROM, POS, ID, REF and ANN."
            )

        if inputs.variants.entity.id != inputs.bam.entity.id:
            self.error(
                "Sample ids of input annotated variants and input bam file do not match. "
                f"Annotated variants have sample id {inputs.variants.entity.id}, while bam "
                f"file has sample id {inputs.bam.entity.id}."
            )

        if inputs.variants.output.species != inputs.reference.output.species:
            self.error(
                "Species for variants and reference file do not match. "
                f"Variants are from {inputs.variants.output.species}, while reference is "
                f"from {inputs.reference.output.species}."
            )

        if inputs.variants.output.build != inputs.reference.output.build:
            self.error(
                "Genome build for variants and reference file do not match. "
                f"Variants have build {inputs.variants.output.build}, while reference has "
                f"build {inputs.reference.output.build}."
            )
        reference = pd.read_csv(inputs.reference.output.tsv.path, sep="\t", header=0)
        check_reference(reference=reference, error=self.error)

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

        variants_table = prepare_variants_table(
            variants_table=variants_table,
            vcf_fields=inputs.vcf_fields,
            ann_fields=inputs.ann_fields,
            gt_fields=inputs.advanced.gf_fields,
            warning=self.warning,
        )
        reference = ann_field_to_df(
            variants_table=reference,
            ann_fields=[
                "Gene_Name",
                "HGVS.p",
            ],
            vcf_fields=["CHROM", "POS", "ID", "REF"],
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
        for field in inputs.vcf_fields:
            output_table[field] = []

        genes = []
        output_table, genes = get_output_table(
            mutations=mutations,
            reference=reference,
            variants_table=variants_table,
            output_table=output_table,
            genes=genes,
            error=self.error,
        )

        get_depth(variants_table=output_table, bam=bam)
        output_table.reset_index(inplace=True, drop=True)
        output_table.to_csv("mutations.tsv", sep="\t", index=False)

        outputs.tsv = "mutations.tsv"
        outputs.genes = genes
        outputs.species = inputs.variants.output.species
        outputs.build = inputs.variants.output.build

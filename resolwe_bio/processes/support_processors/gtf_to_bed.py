"""GTF to BED process."""

import pandas as pd

from resolwe.process import (
    BooleanField,
    DataField,
    FileField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class GTFtoBED(Process):
    """GTF to BED conversion for predefined genes and feature types.

    Note that this process only works with ENSEMBL annotations.
    """

    slug = "gtf-to-bed"
    name = "GTF to BED"
    process_type = "data:bed"
    version = "1.2.0"
    category = "Other"
    data_name = "{{ geneset|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {"cores": 1, "memory": 16384, "storage": 10},
    }

    class Input:
        """Input fields."""

        annotation = DataField(
            "annotation:gtf",
            label="Annotation file (GTF)",
            description="The input GTF file to convert to BED file.",
        )

        feature_type = StringField(
            label="Feature type",
            description="Feature type to extract from the GTF file.",
            default="gene",
            choices=[
                ("gene", "Genes"),
                ("transcript", "Transcripts"),
                ("exon", "Exons"),
            ],
        )

        annotation_source = ListField(
            StringField(),
            label="Annotation source",
            description="Annotation source to extract from the GTF file.",
            default=[
                "ensembl",
                "ensembl_havana",
                "ensembl_havana_tagene",
                "ensembl_tagene",
                "havana",
                "havana_tagene",
                "insdc",
                "mirbase",
            ],
        )

        annotation_field = StringField(
            label="Annotation field",
            description="Annotation field to use as name column in the output BED file. "
            "Note that this field is dependent on the selected feature type. "
            "For gene feature type only gene ID or symbol/name can be extracted. "
            "For transcript feature type only gene ID, gene symbol/name or transcript ID can be extracted. "
            "For exon feature type all the possible fields can be extracted. "
            "Gene name/ID and feature ID will extract gene name/ID and exon/transcript ID, based on the selection of Feature type field. "
            "The IDs will be delimited by '|'. For example ENSG00000279493|ENST00000624081."
            "If exons are selected as a feature type, the name column will include both transcript and exon IDs. "
            "Using None will not include the name column in the BED file or will denote value '.' if option Output strand is selected.",
            default="None",
            choices=[
                ("gene_name", "Gene name"),
                ("gene_id", "Gene ID"),
                ("transcript_id", "Transcript ID"),
                ("exon_id", "Exon ID"),
                ("gene_name_feature_id", "Gene name and feature ID"),
                ("gene_id_feature_id", "Gene ID and feature ID"),
                ("None", "None"),
            ],
        )

        geneset = DataField(
            "geneset",
            label="Gene set",
            description="Gene set to use for filtering.",
            required=True,
        )

        canonical_transcripts = DataField(
            "geneset",
            label="Canonical transcripts",
            description="Canonical transcripts to use for filtering. Only used for transcript and exon feature types.",
            required=False,
            disabled="feature_type == 'gene'",
        )

        output_strand = BooleanField(
            label="Output strand",
            description="Include strand information in the output BED file.",
            default=False,
        )

    class Output:
        """Output fields."""

        bed = FileField(label="BED file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""

        annotation_field = inputs.annotation_field
        feature_type = inputs.feature_type

        if inputs.annotation.output.source != "ENSEMBL":
            self.error("Only ENSEMBL annotations are supported.")

        if feature_type == "genes" and not annotation_field in [
            "gene_id",
            "gene_name",
        ]:
            self.error(
                f"For gene feature type only gene ID or symbol/name can be extracted. You chose {annotation_field}."
            )
        if feature_type == "transcripts" and not annotation_field in [
            "gene_id",
            "gene_name",
            "transcript_id",
            "gene_name_feature_id",
            "gene_id_feature_id",
        ]:
            self.error(
                f"For transcript feature type only gene ID, gene symbol/name or transcript ID can be extracted. You chose {annotation_field}."
            )

        gtf = pd.read_csv(
            inputs.annotation.output.annot.path,
            sep="\t",
            comment="#",
            header=None,
            low_memory=False,
        )

        gtf = pd.DataFrame(
            {
                "chromosome": gtf[0],
                "source": gtf[1],
                "feature_type": gtf[2],
                "start": gtf[3],
                "end": gtf[4],
                "strand": gtf[6],
                "gene_id": gtf[8].str.extract(r'gene_id "(.*?)"', expand=False).values,
                "exon_id": gtf[8].str.extract(r'exon_id "(.*?)"', expand=False).values,
                "transcript_id": gtf[8]
                .str.extract(r'transcript_id "(.*?)"', expand=False)
                .values,
                "gene_name": gtf[8]
                .str.extract(r'gene_name "(.*?)"', expand=False)
                .values,
            },
        )

        gtf = gtf[gtf["source"].isin(inputs.annotation_source)]
        gtf = gtf[gtf["feature_type"] == feature_type]

        if inputs.annotation.output.species != inputs.geneset.output.species:
            self.error(
                "Species of the gene set data object does not match the species of the annotation data object."
            )
        geneset = pd.read_csv(
            inputs.geneset.output.geneset.path,
            delimiter="\t",
            names=["ID"],
        )
        gtf = gtf[gtf["gene_id"].isin(geneset["ID"])]

        if inputs.canonical_transcripts and not feature_type == "gene":
            if (
                inputs.annotation.output.species
                != inputs.canonical_transcripts.output.species
            ):
                self.error(
                    "Canonical transcripts data object species does not match the annotation species."
                )
            canonical_transcripts = pd.read_csv(
                inputs.canonical_transcripts.output.geneset.path,
                sep="\t",
                names=["ID"],
            )
            gtf = gtf[gtf["transcript_id"].isin(canonical_transcripts["ID"])]

        if gtf.empty:
            self.error("No features found with selected inputs.")

        bed = pd.DataFrame(
            {
                "chromosome": gtf["chromosome"],
                "start": gtf["start"] - 1,
                "end": gtf["end"],
            }
        )

        if annotation_field in ["gene_name_feature_id", "gene_id_feature_id"]:
            if feature_type == "exon":
                feature_id = gtf["transcript_id"] + "|" + gtf["exon_id"]
            elif feature_type == "transcript":
                feature_id = gtf["transcript_id"]

            if "gene_name" in annotation_field:
                bed["name"] = gtf["gene_name"] + "|" + feature_id
            elif "gene_id" in annotation_field:
                bed["name"] = gtf["gene_id"] + "|" + feature_id
        else:
            bed["name"] = gtf[annotation_field] if annotation_field != "None" else "."

        if inputs.output_strand:
            bed["score"] = "."
            bed["strand"] = gtf["strand"]
        else:
            if annotation_field == "None":
                bed = bed.drop(columns=["name"])

        sorted_bed = bed.sort_values(by=["chromosome", "start", "end"])
        sorted_bed.drop_duplicates(inplace=True)

        bed_fn = "sorted_bed.bed"

        sorted_bed.to_csv(bed_fn, sep="\t", header=False, index=False)

        outputs.bed = bed_fn
        outputs.species = inputs.annotation.output.species
        outputs.build = inputs.annotation.output.build

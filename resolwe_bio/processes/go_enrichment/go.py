"""Gene Ontology Enrichment Analysis."""

import csv
import tempfile
from collections import defaultdict

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FloatField,
    IntegerField,
    JsonField,
    ListField,
    Persistence,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


class GOEnrichmentAnalysis(ProcessBio):
    """Identify significantly enriched Gene Ontology terms for given genes."""

    slug = "goenrichment"
    name = "GO Enrichment analysis"
    process_type = "data:goea"
    version = "3.6.3"
    category = "Enrichment and Clustering"
    data_name = 'GO Enrichment analysis for {{genes|join(", ")|default("?")}}'
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {"cores": 1, "memory": 4096, "storage": 10},
    }

    class Input:
        """Input fields to process GOEnrichmentAnalysis."""

        ontology = DataField("ontology:obo", label="Gene Ontology")
        gaf = DataField("gaf", label="GO annotation file (GAF v2.0)")
        genes = ListField(
            StringField(), label="List of genes", placeholder="new gene id"
        )
        source = StringField(label="Gene ID database of selected genes")
        species = StringField(
            label="Species",
            allow_custom_choice=True,
            description="Specify species name. This field is required if gene subset is set.",
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
            ],
        )
        pval_threshold = FloatField(
            label="P-value threshold",
            default=0.1,
        )
        min_genes = IntegerField(
            label="Minimum number of genes",
            default=1,
            description="Minimum number of genes on a GO term.",
        )

    class Output:
        """Output fields to process GOEnrichmentAnalysis."""

        terms = JsonField(label="Enriched terms")
        ids = FileField(label="Mapped ids", required=False)
        source = StringField(label="Source")
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run analysis."""

        if inputs.species != inputs.gaf.output.species:
            self.warning(
                "Selected genes Species must be the same as the Species field of the GAF file."
            )
            self.error(
                f"Selected genes are from {inputs.species}, "
                f"while GAF file has defined {inputs.gaf.output.species} under Species field."
            )

        org_features = self.feature.filter(
            source=inputs.source, species=inputs.species, feature_id__in=inputs.genes
        )

        if len(org_features) == 0:
            self.error("No genes were fetched from the knowledge base.")

        if inputs.source == inputs.gaf.output.source:
            target_ids = inputs.genes
        else:
            mapping_res = self.mapping.filter(
                source_db=inputs.source,
                source_species=inputs.species,
                target_db=inputs.gaf.output.source,
                target_species=inputs.gaf.output.species,
                source_id__in=inputs.genes,
            )

            if len(mapping_res) == 0:
                self.error("Failed to map features.")

            ids = defaultdict(list)
            target_ids = []
            for m in mapping_res:
                if m.source_id in inputs.genes:
                    target_ids.append(m.target_id)
                    if m.source_id in ids:
                        self.warning(f"Mapping {m.source_id} returned multiple times.")
                    ids[m.source_id].append(m.target_id)

            if len(inputs.genes) > len(ids):
                self.warning("Not all features could be mapped.")

            if len(target_ids) > 0:
                with open("mapped_ids.txt", "w") as f:
                    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                    writer.writerow([inputs.source, inputs.gaf.output.source])
                    for key, value in ids.items():
                        for v in value:
                            writer.writerow([key, v])
                outputs.ids = "mapped_ids.txt"

        with tempfile.NamedTemporaryFile() as input_genes:
            input_genes.write(" ".join(target_ids).encode("UTF-8"))
            input_genes.flush()

            args = [
                str(inputs.pval_threshold),
                str(inputs.min_genes),
                inputs.ontology.output.obo_obj.path,
                inputs.gaf.output.gaf_obj.path,
                input_genes.name,
            ]

            (Cmd["processor"][args] > "terms.json")()

        outputs.source = inputs.gaf.output.source
        outputs.species = inputs.gaf.output.species
        outputs.terms = "terms.json"

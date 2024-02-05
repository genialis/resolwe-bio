"""Upload ML-ready expression matrix."""

from pathlib import Path

import pandas as pd

from resolwe.process import DataField, FileField, Process, SchedulingClass, StringField
from resolwe.process.models import Entity


class UploadMLExpression(Process):
    """Upload ML-ready expression matrix."""

    slug = "upload-ml-expression"
    name = "ML-ready expression"
    process_type = "data:ml:table:expressions"
    version = "1.0.2"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
            "storage": 10,
        },
    }
    data_name = "{{ reference_space|name }}"

    class Input:
        """Inputs."""

        exp = FileField(
            label="Transformed expressions",
            description="A TAB separated file containing transformed "
            "expression values with sample IDs for index (first column "
            "with label sample_id) and ENSEMBL IDs (recommended but not "
            "required) for the column names.",
        )
        source = StringField(
            label="Feature source",
            allow_custom_choice=True,
            choices=[
                ("AFFY", "AFFY"),
                ("DICTYBASE", "DICTYBASE"),
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )
        species = StringField(
            label="Species",
            description="Species latin name.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
            ],
        )
        reference_space = DataField(
            "ml:space", label="Reference space of ML-ready data"
        )

    class Output:
        """Outputs."""

        exp = FileField(label="Transformed expressions")
        source = StringField(label="Feature source")
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run the analysis."""

        # Parse exp file:
        exp_path = inputs.exp.import_file(imported_format="extracted")
        try:
            df = pd.read_csv(
                exp_path,
                sep="\t",
                index_col=0,
                float_precision="round_trip",
            )
        except Exception as err:
            self.error(f"It was not possible to read the provided table. {err}")

        # Sort columns in alphabetical order
        df.sort_index(axis=1, inplace=True)

        # Ensure that index contains actual ids from expressions
        sample_ids_df = df.index.astype(int).tolist()
        samples = Entity.filter(id__in=sample_ids_df)
        sample_ids_server = [s.id for s in samples]
        missing = set(sample_ids_df) - set(sample_ids_server)
        if missing:
            missing_str = ", ".join(map(str, list(missing)[:5])) + (
                "..." if len(missing) > 5 else ""
            )
            self.error(
                f"There are {len(missing)} samples in uploaded matrix that "
                f"are inaccessible or missing on server: {missing_str}"
            )

        # Ensure that sources of this object and reference space are equal:
        if inputs.source != inputs.reference_space.output.source:
            self.error(
                "Source of expression matrix must match the source of reference space."
            )
        # Ensure that species of this object and reference space are equal:
        if inputs.species != inputs.reference_space.output.species:
            self.error(
                "Species of expression matrix must match the species of reference space."
            )

        # Ensure that expression has the same features as the reference space
        features_ref = inputs.reference_space.output.features.json["features"]
        # Both should be sorted by now, so one can directly compare:
        if features_ref != df.columns.tolist():
            self.error(
                "Uploaded expression matrix does not have the same features as the reference space."
            )

        # Save
        final_exp_name = Path(exp_path).name
        df.to_csv(final_exp_name, sep="\t")

        # Set outputs
        outputs.exp = final_exp_name
        outputs.source = inputs.source
        outputs.species = inputs.species

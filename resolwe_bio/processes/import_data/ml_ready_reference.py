"""Upload reference space for ML ready data."""

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from resolwe.process import FileField, JsonField, SchedulingClass, StringField
from resolwe.process.models import Collection, DescriptorSchema, Entity

from resolwe_bio.process.runtime import ProcessBio


class ReferenceSpace(ProcessBio):
    """
    Define the reference space of ML-ready data sets.

    Use a descriptive name that uniquely defines the reference space, e.g.:

        - ComBat space of all TCGA v001
        - Quantile transform to uniform distribution v042.
    """

    slug = "reference-space"
    name = "Reference space"
    process_type = "data:ml:space"
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
    data_name = "{{ name }}"

    class Input:
        """Inputs."""

        name = StringField(label="Reference space name")
        description = StringField(label="Reference space description")
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
        training_data = FileField(
            label="Traning data",
            description="A TAB separated file containing expression values "
            "used to create the preprocessor. The data should have Sample "
            "object ID for index (first column with label sample_id) and "
            "Ensembl gene IDs (recommended but not required) for the "
            "column names.",
        )
        preprocessor = FileField(
            label="Pickled preprocessor",
            description="Serialized (pickled) preprocessor used to transform "
            "data to the reference space.",
            required=False,
        )

    class Output:
        """Outputs."""

        features = JsonField(label="List of features")
        source = StringField(label="Feature ID source")
        species = StringField(label="Species")
        training_data = FileField(label="Traning data")
        preprocessor = FileField(label="Pickled preprocessor")

    def run(self, inputs, outputs):
        """Run the analysis."""

        # Ensure the object is uploaded to the general collection
        reference_spaces_collection_slug = "reference-spaces"
        collections = Collection.filter(data__id=self.data.id)
        if not collections:
            self.warning("Reference space table was not uploaded to a Collection.")
        elif collections[0].slug != reference_spaces_collection_slug:
            self.warning(
                f"Reference space table was not uploaded to Collection {reference_spaces_collection_slug}."
            )

        # Parse file:
        training_data_file = inputs.training_data.import_file(
            imported_format="extracted"
        )
        try:
            df = pd.read_csv(
                training_data_file,
                sep="\t",
                float_precision="round_trip",
                index_col=0,
            )
        except Exception as err:
            self.error(f"It was not possible to read the provided table. {err}")

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
                f"There are {len(missing)} samples in uploaded table that "
                f"are missing on server: {missing_str}"
            )

        features_df = sorted(df.columns.tolist())
        # If source is in KB, ensure that columns contain features from provided source
        if not self.feature.exists(source=inputs.source):
            self.warning(f"There are no features with source={inputs.source} in KB.")
        else:
            features_kb = self.feature.filter(
                feature_id__in=features_df,
                species=inputs.species,
            )
            source_counter = Counter([f.source for f in features_kb])
            sources = list(source_counter.keys())
            if len(sources) > 1:
                sources = ", ".join(sources)
                self.error(
                    f"Features in training data are from different sources: {sources}"
                )
            elif sources[0] != inputs.source:
                self.error(
                    "Features in training data are not from the source specified in inputs."
                )
            elif source_counter[inputs.source] != len(features_df):
                self.error("Some features in training data are not in KB")

        # Sort columns and index and save to file
        df = df.sort_index(axis=0).sort_index(axis=1)
        training_data_file = Path(training_data_file).name
        df.to_csv(training_data_file, sep="\t")

        # Extract feature/genes from expressions
        features_json_name = "features.json"
        with open(features_json_name, "wt") as handle:
            json.dump({"features": features_df}, handle)

        # Get preprocessor pickle file
        preprocessor_file = inputs.preprocessor.import_file(imported_format="extracted")

        # Set the descriptor schema to geneset and attach description
        ds = DescriptorSchema.get_latest(slug="geneset")
        self.data.descriptor_schema = ds.id
        self.data.descriptor = {"description": inputs.description}

        # Set outputs
        outputs.source = inputs.source
        outputs.species = inputs.species
        outputs.features = features_json_name
        outputs.preprocessor = preprocessor_file
        outputs.training_data = training_data_file

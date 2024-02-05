"""Upload metadata table."""

from pathlib import Path

import pandas as pd

from resolwe.process import FileField, IntegerField, Process, SchedulingClass
from resolwe.process.models import Collection, Entity

SAMPLE_COLUMNS = {
    "Sample ID": "id",
    "Sample slug": "slug",
    "Sample name": "name",
}


def lower_suffix(path, info):
    """Change suffix of a file to lowercase."""
    if not path.suffix.islower():
        new_path = path.with_suffix(path.suffix.lower())
        path.replace(new_path)
        info("File extension of the table was replaced with a lower case version.")
        return new_path
    else:
        return path


def read_tabular_data(path, sample_columns, error):
    """Convert the uploaded file to Pandas data frame."""
    extensions = [".csv", ".tab", ".tsv", ".xlsx", ".xls"]
    if path.suffix not in extensions:
        error(
            "Unsupported file name extension. Supported extensions "
            f"are {', '.join(extensions)}."
        )
    try:
        if path.suffix == ".xls":
            df = pd.read_excel(path, engine="xlrd")
        elif path.suffix == ".xlsx":
            df = pd.read_excel(path, engine="openpyxl")
        elif any(path.suffix == ext for ext in [".tab", ".tsv"]):
            df = pd.read_csv(path, sep="\t")
        elif path.suffix == ".csv":
            df = pd.read_csv(path)
        else:
            df = pd.DataFrame()
    except Exception as err:
        error(f"It was not possible to read the provided data table. {err}")

    if len(df.columns.intersection(sample_columns)) != 1:
        error(
            f"The uploaded metadata table needs to contain "
            f"exactly one of the following columns: "
            f"{sorted(sample_columns.keys())}."
        )

    if len(df) < 1:
        error("The uploaded table contains no samples.")

    return df


class UploadMetadataUnique(Process):
    """Upload metadata file where each row corresponds to a single sample.

    The uploaded metadata table represents one-to-one (1:1) relation to
    samples in the working collection. Metadata table must contain a column
    with one of the following headers: "Sample ID", "Sample name" or "Sample slug".
    """

    slug = "upload-metadata-unique"
    name = "Metadata table (one-to-one)"
    process_type = "data:metadata:unique"
    version = "1.1.1"
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
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input field to process UploadMetadataUnique."""

        src = FileField(
            label="Table with metadata",
            description="The metadata table should use one of the following "
            "extensions: .csv, .tab, .tsv, .xlsx, .xls",
        )

    class Output:
        """Output field of the process UploadMetadataUnique."""

        table = FileField(label="Uploaded table")
        n_samples = IntegerField(label="Number of samples")

    def run(self, inputs, outputs):
        """Run the analysis."""

        collections = Collection.filter(data__id=self.data.id)
        if not collections:
            self.error(
                "Metadata table was not uploaded to a Collection. "
                "Matching of metadata entries to Sample objects is not possible."
            )

        samples = Entity.filter(collection_id=collections[0].id)

        path = Path(inputs.src.import_file(imported_format="extracted"))

        # change the file suffix if it is either upper or mixed case
        path = lower_suffix(path, info=self.info)

        df_data = read_tabular_data(path, SAMPLE_COLUMNS, error=self.error)

        sample_header = df_data.columns.intersection(SAMPLE_COLUMNS)[0]

        col_samples = {
            getattr(sample, SAMPLE_COLUMNS[sample_header]) for sample in samples
        }

        df_samples = df_data[sample_header]
        intersection = col_samples.intersection(df_samples.values)

        if not intersection:
            self.warning(
                "None of the samples listed in the uploaded Sample metadata table "
                "match the Samples in the working Collection."
            )

        dup_samples = df_samples[df_samples.duplicated()]
        if not dup_samples.empty:
            self.error(
                f"Duplicated metadata entries {dup_samples.tolist()} were found. "
                f"Please use the metadata upload process that "
                f"allows for one-to-many relations instead."
            )

        outputs.table = str(path)
        outputs.n_samples = len(df_samples.unique())


class UploadMetadata(Process):
    """Upload metadata file where more than one row can match to a single sample.

    The uploaded metadata table represents one-to-many (1:n) relation to
    samples in the working collection. Metadata table must contain a column
    with one of the following headers: "Sample ID", "Sample name" or "Sample slug".
    """

    slug = "upload-metadata"
    name = "Metadata table"
    process_type = "data:metadata"
    version = "1.1.1"
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
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input field to process UploadMetadata."""

        src = FileField(
            label="Table with metadata",
            description="The metadata table should use one of the following "
            "extensions: .csv, .tab, .tsv, .xlsx, .xls",
        )

    class Output:
        """Output field of the process UploadMetadata."""

        table = FileField(label="Uploaded table")
        n_samples = IntegerField(label="Number of samples")

    def run(self, inputs, outputs):
        """Run the analysis."""

        collections = Collection.filter(data__id=self.data.id)
        if not collections:
            self.error(
                "Metadata table was not uploaded to a Collection. "
                "Matching of metadata entries to Sample objects is not possible."
            )

        samples = Entity.filter(collection_id=collections[0].id)

        path = Path(inputs.src.import_file(imported_format="extracted"))

        # change the file suffix if it is either upper or mixed case
        path = lower_suffix(path, info=self.info)

        df_data = read_tabular_data(path, SAMPLE_COLUMNS, error=self.error)

        sample_header = df_data.columns.intersection(SAMPLE_COLUMNS)[0]

        col_samples = {
            getattr(sample, SAMPLE_COLUMNS[sample_header]) for sample in samples
        }

        df_samples = df_data[sample_header]
        intersection = col_samples.intersection(df_samples.values)

        if not intersection:
            self.warning(
                "None of the samples listed in the uploaded Sample metadata table "
                "match the Samples in the working Collection."
            )

        outputs.table = str(path)
        outputs.n_samples = len(df_samples.unique())

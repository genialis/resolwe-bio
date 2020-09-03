"""Upload metadata table in Orange format."""
from pathlib import Path

import Orange

from resolwe.process import (
    FileField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


def format_missing(missing):
    """Format number of missing values as percentage."""
    if missing:
        return f"({missing:.1%} missing values)"
    else:
        return "(no missing values)"


def get_features_description(data):
    """Get formatted features description."""
    missing_features = format_missing(
        data.has_missing_attribute() and data.get_nan_frequency_attribute()
    )
    return f"{len(data.domain.attributes)} {missing_features}"


def get_target_description(data):
    """Get formatted target class description."""
    missing_in_class = format_missing(
        data.has_missing_class() and data.get_nan_frequency_attribute()
    )

    if data.domain.has_continuous_class:
        target_description = f"Regression; numerical class {missing_in_class}"
    elif data.domain.has_discrete_class:
        target_description = (
            "Classification; categorical class "
            f"with {len(data.domain.class_var.values)} values {missing_in_class}"
        )
    elif data.domain.class_vars:
        target_description = (
            "Multi-target; "
            f"{len(data.domain.class_vars)} target variables "
            f"{missing_in_class}"
        )
    else:
        target_description = "Not defined"

    return target_description


def change_suffix(path):
    """Change suffix of a file to lowercase."""
    new_path = path.with_suffix(path.suffix.lower())
    path.replace(new_path)
    return new_path


class UploadOrangeMetadata(Process):
    """Upload metadata table in Orange format.

    Orange can read files in native tab-delimited format, or can load
    data from any of the major standard spreadsheet file types, like CSV
    and Excel. Native format starts with the names of attributes with
    prefixes that define attribute type (continuous, discrete, time,
    string) and role (class, meta, ignore, instance weights). Prefixes
    are separated from the attribute name with a hash sign (“#”).

    Legacy format with three header rows is also supported. The first
    row lists feature (column) names. The second header row gives the
    attribute type and the third header line contains role information
    to identify dependent features (class), irrelevant features (ignore)
    or meta features (meta).

    For more information see Orange
    [documentation](https://orange-visual-programming.readthedocs.io/loading-your-data/index.html#header-with-attribute-type-information).

    An example of native tab-delimited format can be downloaded
    [here](http://file.biolab.si/datasets/sample-head.xlsx).

    """

    slug = "upload-orange-metadata"
    name = "Metadata table for Orange"
    process_type = "data:metadata:orange"
    version = "1.0.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/orange:1.0.0"}},
        "resources": {"cores": 1, "memory": 8192},
    }
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input field to process UploadOrangeMetadata."""

        src = FileField(
            label="Table with metadata",
            description="The table should be in Orange format and use "
            "one of the following extensions: .csv, .tab, .tsv, .xlsx, "
            ".xls",
        )

    class Output:
        """Output field of the process UploadOrangeMetadata."""

        table = FileField(label="Uploaded table")
        n_samples = IntegerField(label="Number of samples")
        features = StringField(label="Number of features")
        target = StringField(label="Target class description")
        n_metas = IntegerField(label="Number of meta attributes")

    def run(self, inputs, outputs):
        """Run the analysis."""

        table_path = inputs.src.import_file(imported_format="extracted")

        extensions = [".csv", ".tab", ".tsv", ".xlsx", ".xls"]
        path = Path(table_path)

        if path.suffix in [e.upper() for e in extensions]:
            path = change_suffix(path)
            self.info(
                "File extension of the table was replaced with a lower case version."
            )

        if path.suffix not in extensions:
            self.error(
                "Unsupported file name extension. Supported extensions "
                f"are {', '.join(extensions)}."
            )

        try:
            data = Orange.data.Table.from_file(str(path))
        except Exception as err:
            outputs.table = str(path)
            self.error(f"Orange is unable to read the provided data table. {err}")

        n_samples = len(data)
        if n_samples < 1:
            self.error("The uploaded table contains no samples.")

        outputs.table = str(path)
        outputs.n_samples = n_samples
        outputs.features = get_features_description(data)
        outputs.target = get_target_description(data)
        outputs.n_metas = len(data.domain.metas)

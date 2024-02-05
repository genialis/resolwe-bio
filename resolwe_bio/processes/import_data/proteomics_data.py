"""Upload proteomics data."""

import re
from pathlib import Path

import pandas as pd

from resolwe.process import FileField, Process, SchedulingClass, StringField


def change_suffix(path):
    """Change suffix of a file to lowercase."""
    new_path = path.with_suffix(path.suffix.lower())
    path.replace(new_path)
    return new_path


def prepare_filename(fname):
    """Return a sanitized string that can be used as a file name."""
    return re.sub(r"(?u)[^-\w.]", "", str(fname).strip().replace(" ", "_"))


class UploadProteomicsData(Process):
    """Upload a mass spectrometry proteomics sample data file.

    The input 5-column tab-delimited file with the .txt suffix is
    expected to contain a header line with the following meta-data
    column names: "Uniprot ID", "Gene symbol", "Protein name" and
    "Number of peptides". The fifth column contains the sample data.
    """

    slug = "upload-proteomics-sample"
    name = "Upload proteomics sample"
    process_type = "data:proteomics:massspectrometry"
    version = "1.2.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {"cores": 1, "memory": 2048},
    }
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
    }
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input field to process UploadProteomicsData."""

        src = FileField(label="Table containing mass spectrometry data (.txt)")
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu "
            "or write a custom species name in the species field.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
            ],
        )
        source = StringField(
            label="Protein ID database source",
            default="UniProtKB",
            choices=[("UniProtKB", "UniProtKB")],
        )

    class Output:
        """Output field of the process UploadProteomicsData."""

        table = FileField(label="Uploaded table")
        species = StringField(label="Species")
        source = StringField(label="Source")

    def run(self, inputs, outputs):
        """Run the analysis."""

        table_path = inputs.src.import_file(imported_format="extracted")

        extensions = [".txt"]
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

        required_columns = [
            "Uniprot ID",
            "Gene symbol",
            "Protein name",
            "Number of peptides",
        ]
        table = pd.read_csv(path, sep="\t")
        header = list(table)
        if not set(required_columns).issubset(header):
            self.error(
                f"The input file must contain all of the required columns: {required_columns}."
            )
        if not len(header) == 5:
            self.error(
                f"The input file must contain the required metadata columns: {required_columns} "
                f"and exactly one sample-data column. The provided input file contains columns: "
                f"{header}."
            )

        outputs.table = str(path)
        outputs.species = inputs.species
        outputs.source = inputs.source


class UploadProteomicsDataSet(Process):
    """Upload a mass spectrometry proteomics sample set file.

    The input multi-sample tab-delimited file with the .txt suffix is
    expected to contain a header line with the following meta-data
    column names: "Uniprot ID", "Gene symbol", "Protein name" and
    "Number of peptides". Each additional column in the input file
    should contain data for a single sample.
    """

    slug = "upload-proteomics-sample-set"
    name = "Upload proteomics sample set"
    process_type = "data:proteomics:sampleset"
    version = "1.2.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {"cores": 1, "memory": 2048},
    }
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input field to process UploadProteomicsDataSet."""

        src = FileField(label="Table containing mass spectrometry data (.txt)")
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu "
            "or write a custom species name in the species field.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
            ],
        )
        source = StringField(
            label="Protein ID database source",
            default="UniProtKB",
            choices=[("UniProtKB", "UniProtKB")],
        )

    class Output:
        """Output field of the process UploadProteomicsDataSet."""

        table = FileField(label="Uploaded table")
        species = StringField(label="Species")
        source = StringField(label="Source")

    def run(self, inputs, outputs):
        """Run the analysis."""

        table_path = inputs.src.import_file(imported_format="extracted")

        extensions = [".txt"]
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

        outputs.table = str(path)
        outputs.species = inputs.species
        outputs.source = inputs.source

        # spawn individual samples from the input sample set file
        required_columns = [
            "Gene symbol",
            "Protein name",
            "Number of peptides",
        ]
        sample_set = pd.read_csv(path, sep="\t", index_col="Uniprot ID")
        header = list(sample_set)

        if not set(required_columns).issubset(header):
            self.error(
                f"The input file must contain all of the required columns: {required_columns}."
            )

        for sample_column in header:
            if sample_column not in required_columns:
                sample_data = sample_set[required_columns + [sample_column]]
                sample_data_name = prepare_filename(sample_column) + ".txt"
                sample_data.to_csv(sample_data_name, index_label="Uniprot ID", sep="\t")

                # spawn a new sample as an individual object
                process_inputs = {
                    "src": sample_data_name,
                    "species": inputs.species,
                    "source": inputs.source,
                }
                self.run_process("upload-proteomics-sample", process_inputs)

"""Upload microarray expression data."""

import gzip
import io
import json
from pathlib import Path

import pandas as pd

from resolwe.process import (
    DataField,
    FileField,
    JsonField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


def gzopen(fname):
    """Open Gzip files using io.BufferedReader."""
    return io.TextIOWrapper(io.BufferedReader(gzip.open(fname)))


def isfloat(value):
    """Check if value is float."""
    try:
        float(value)
        return True
    except ValueError:
        return False


def expression_to_json(infile, outfile):
    """Re-format expression file to json."""
    with gzopen(infile) as f:
        # Split lines by tabs
        # Ignore lines without a number in second column
        # Build a dictionary of gene-expression pairs
        exp = {
            "genes": {
                gene_exp[0]: float(gene_exp[1])
                for gene_exp in (l.split("\t") for l in f)
                if len(gene_exp) == 2 and isfloat(gene_exp[1])
            }
        }

        with open(outfile, "w") as f:
            json.dump(exp, f)


def rename_cols(infile, outfile):
    """Rename columns in expression file."""
    exp = pd.read_csv(
        infile,
        compression="gzip",
        sep="\t",
        skip_blank_lines=True,
        usecols=["Gene", "Expression"],
        index_col="Gene",
        dtype={
            "Gene": str,
            "Expression": float,
        },
    )
    exp.squeeze()
    return exp.to_csv(
        outfile,
        index_label="FEATURE_ID",
        header=["log2 normalized intensity signal"],
        sep="\t",
    )


def prepare_expression_set(exp_file, feature_dict, outfile_name):
    """Prepare expression set output data."""
    exp = pd.read_csv(exp_file, sep="\t", float_precision="round_trip")
    exp["FEATURE_ID"] = exp["FEATURE_ID"].astype("str")
    exp["GENE_SYMBOL"] = exp["FEATURE_ID"].map(feature_dict)
    input_features = exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    # Reorder columns
    columns = ["FEATURE_ID", "GENE_SYMBOL", "log2 normalized intensity signal"]
    exp_set = exp[columns]
    # Replace NaN values with empty string
    exp_set.fillna("", inplace=True)

    # Write to file
    exp_set.to_csv(
        outfile_name + ".txt.gz",
        header=True,
        index=False,
        sep="\t",
        compression="gzip",
    )

    # Write to JSON
    df_dict = exp_set.set_index("FEATURE_ID").to_dict(orient="index")
    with open(outfile_name + ".json", "w") as f:
        json.dump({"genes": df_dict}, f, allow_nan=False)


class ImportMicroarrayExpression(Process):
    """Import unmapped microarray expression data."""

    slug = "upload-microarray-expression"
    name = "Upload microarray expression (unmapped)"
    process_type = "data:microarray:normalized"
    version = "1.1.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 4096,
            "network": True,
        },
    }
    data_name = '{{ exp.file|default("?") }}'

    class Input:
        """Input fields to process ImportMicroarrayExpression."""

        exp = FileField(
            label="Normalized expression",
            description="Normalized expression file with the original probe IDs. Supported file extensions are "
            ".tab.*, .tsv.*, .txt.*",
        )
        exp_type = StringField(
            label="Normalization type",
        )
        platform = StringField(
            label="Microarray platform name",
        )
        platform_id = StringField(
            label="GEO platform ID",
            description="Platform ID according to the GEO database. This can be used in following steps to "
            "automatically map probe IDs to genes.",
            required=False,
        )
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu or write a custom species name in the species "
            "field",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Macaca mulatta", "Macaca mulatta"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
            ],
        )

    class Output:
        """Output fields to process ImportMicroarrayExpression."""

        exp = FileField(label="Uploaded normalized expression")
        exp_type = StringField(label="Normalization type")
        platform = StringField(label="Microarray platform type")
        platform_id = StringField(label="GEO platform ID", required=False)
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run the analysis."""
        exp = inputs.exp.import_file(imported_format="compressed")

        supported_extensions = (".tab", ".tsv", ".txt")
        if not Path(exp).stem.endswith(supported_extensions):
            self.error(
                f"The imported file has unsupported file name extension. "
                f"The supported extensions are {supported_extensions}."
            )

        if inputs.platform_id:
            outputs.platform_id = inputs.platform_id

        outputs.exp = exp
        outputs.exp_type = inputs.exp_type
        outputs.platform = inputs.platform
        outputs.species = inputs.species


class MicroarrayExpression(ProcessBio):
    """Upload normalized and mapped microarray expression data."""

    slug = "mapped-microarray-expression"
    name = "Mapped microarray expression"
    process_type = "data:expression:microarray"
    version = "1.3.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 4096,
            "network": True,
        },
    }
    data_name = "{{ exp_unmapped|name|default('?') }}"

    class Input:
        """Input fields to process MicroarrayExpression."""

        exp_unmapped = DataField(
            "microarray:normalized",
            label="Unmapped normalized expressions",
            description="Unmapped normalized expression with the original probe IDs.",
        )
        exp = FileField(
            label="Normalized and mapped expressions file",
            description="Files should have two columns one with GeneIDs and the other one with expression values."
            "Expected column names are 'Gene' and 'Expression'.Supported file extensions are .tab.*, .tsv.*, .txt.*",
        )
        source = StringField(
            label="Gene ID source",
            allow_custom_choice=True,
            choices=[
                ("AFFY", "AFFY"),
                ("DICTYBASE", "DICTYBASE"),
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )
        build = StringField(
            label="Genome build",
        )
        probe_mapping = StringField(
            label="Probe to transcript mapping used",
        )

    class Output:
        """Output fields."""

        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        platform = StringField(label="Microarray platform type")
        platform_id = StringField(label="GEO platform ID", required=False)
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")
        probe_mapping = StringField(label="Probe to transcript mapping used")

    def run(self, inputs, outputs):
        """Run the analysis."""
        exp = inputs.exp.import_file(imported_format="compressed")
        exp_stem = Path(exp).stem

        supported_extensions = (".tab", ".tsv", ".txt")
        if not exp_stem.endswith(supported_extensions):
            self.error(
                "The imported file has unsupported file name extension. "
                f"The supported extensions are {supported_extensions}."
            )

        name = exp_stem[:-4]

        expression_to_json(exp, "json.txt")

        # Rename columns of the expression file
        exp_renamed = f"{exp}_renamed"
        rename_cols(infile=exp, outfile=exp_renamed)

        # Prepare the expression set outputs
        feature_ids = pd.read_csv(
            exp_renamed, sep="\t", index_col="FEATURE_ID"
        ).index.tolist()

        feature_filters = {
            "source": inputs.source,
            "species": inputs.exp_unmapped.output.species,
            "feature_id__in": feature_ids,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        prepare_expression_set(
            exp_file=exp_renamed,
            feature_dict=feature_ids_to_names,
            outfile_name=f"{name}_expressions",
        )

        if inputs.exp_unmapped.output.platform_id:
            outputs.platform_id = inputs.exp_unmapped.output.platform_id

        outputs.exp = exp
        outputs.exp_json = "json.txt"
        outputs.exp_type = inputs.exp_unmapped.output.exp_type
        outputs.platform = inputs.exp_unmapped.output.platform
        outputs.exp_set = name + "_expressions.txt.gz"
        outputs.exp_set_json = name + "_expressions.json"
        outputs.source = inputs.source
        outputs.species = inputs.exp_unmapped.output.species
        outputs.build = inputs.build
        outputs.feature_type = "gene"
        outputs.probe_mapping = inputs.probe_mapping

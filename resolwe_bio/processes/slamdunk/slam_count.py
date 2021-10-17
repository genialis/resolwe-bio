"""Convert collapsed Slamdunk tcount data to the expressions data type."""
import gzip
import io
import json
import os

import pandas as pd

from resolwe.process import DataField, FileField, JsonField, StringField

from resolwe_bio.process.runtime import ProcessBio


def prepare_expressions(infile, rc_file, tpm_file):
    """Compute TPM values from T>C read counts."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        usecols=["gene_name", "length", "tcReadCount"],
        index_col="gene_name",
        dtype={
            "gene_name": str,
            "length": int,
            "tcReadCount": int,
        },
    )
    exp["rpk"] = exp.apply(lambda x: (x.tcReadCount * 1e3 / x.length), axis=1)
    rpk_sum = exp["rpk"].sum()
    exp["tpm"] = exp.apply(lambda x: (x.rpk / rpk_sum * 1e6), axis=1)
    rc = exp["tcReadCount"].to_csv(
        rc_file,
        index_label="FEATURE_ID",
        header=["RAW_COUNT"],
        compression="gzip",
        sep="\t",
    )
    tpm = exp["tpm"].to_csv(
        tpm_file,
        index_label="Gene",
        header=["Expression"],
        sep="\t",
        compression="gzip",
    )
    return (rc, tpm)


def exp_to_df(exp_file, exp_type):
    """Prepare expression file for gene sets merging."""
    with gzip.open(exp_file) as exp:
        df = pd.read_csv(exp, sep="\t", float_precision="round_trip")
        df.rename(
            index=str,
            columns={"Gene": "FEATURE_ID", "Expression": exp_type},
            inplace=True,
        )
        # Cast FEATURE_ID column to string
        df["FEATURE_ID"] = df["FEATURE_ID"].astype("str")

    return df


def prepare_expression_set(rc, tpm, feature_dict, outfile_name):
    """Prepare expression set output data."""
    rc_exp = pd.read_csv(rc, sep="\t", float_precision="round_trip")
    rc_exp["FEATURE_ID"] = rc_exp["FEATURE_ID"].astype("str")
    tpm_exp = exp_to_df(tpm, "TPM")
    rc_exp["GENE_SYMBOL"] = rc_exp["FEATURE_ID"].map(feature_dict)
    input_features = rc_exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(rc_exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    # Merge with normalized expression values
    exp_set = rc_exp.merge(tpm_exp, on="FEATURE_ID")
    # Reorder columns
    columns = ["FEATURE_ID", "GENE_SYMBOL", "RAW_COUNT", "TPM"]
    exp_set = exp_set[columns]
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


def expression_to_storage(tpm_input, tpm_json):
    """Convert expressions file to JSON format."""

    def isfloat(value):
        """Check if value is float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    with io.TextIOWrapper(io.BufferedReader(gzip.open(tpm_input))) as f:
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

    with open(file=tpm_json, mode="wt") as f:
        json.dump(exp, f)

    return tpm_json


class SlamCount(ProcessBio):
    """Convert gene-level Slamdunk data to the expression data type."""

    slug = "slam-count"
    process_type = "data:expression:slamdunk"
    name = "Slam count"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/slamdunk:2.0.0"},
        },
        "resources": {
            "cores": 1,
            "memory": 4096,
            "network": True,
        },
    }
    category = "Slamdunk"
    entity = {
        "type": "sample",
    }
    data_name = '{{ tcount|sample_name|default("?") }}'
    version = "1.2.0"

    class Input:
        """Input fields for SlamCount."""

        tcount = DataField("alleyoop:collapse", label="Collapsed Slamdunk results")
        source = StringField(
            label="Gene ID source",
            default="ENSEMBL",
            choices=[
                ("ENSEMBL", "ENSEMBL"),
                ("UCSC", "UCSC"),
            ],
        )

    class Output:
        """Output fields to process SlamCount."""

        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        rc = FileField(label="T>C read counts")
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.tcount.output.tcount.path)
        assert basename.endswith(".txt")
        name = basename[:-4]

        rc_file = name + "_rc.txt.gz"
        tpm_file = name + "_tmp.txt.gz"

        prepare_expressions(inputs.tcount.output.tcount.path, rc_file, tpm_file)

        for exp_file in [rc_file, tpm_file]:
            if not os.path.isfile(exp_file):
                self.error(
                    "Failed to parse tcout file. {} file was not created".format(
                        exp_file
                    )
                )

        # Save the abundance estimates to JSON storage
        json_output = "json.txt"
        expression_to_storage(tpm_input=tpm_file, tpm_json=json_output)

        # Prepare the expression set outputs
        feature_ids = pd.read_csv(
            rc_file, sep="\t", compression="gzip", index_col="FEATURE_ID"
        ).index.tolist()

        feature_filters = {
            "source": inputs.source,
            "species": inputs.tcount.output.species,
            "feature_id__in": feature_ids,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        prepare_expression_set(
            rc=rc_file,
            tpm=tpm_file,
            feature_dict=feature_ids_to_names,
            outfile_name=f"{name}_expressions",
        )

        outputs.exp = tpm_file
        outputs.exp_json = json_output
        outputs.exp_type = "TPM"
        outputs.rc = rc_file
        outputs.exp_set = name + "_expressions.txt.gz"
        outputs.exp_set_json = name + "_expressions.json"
        outputs.species = inputs.tcount.output.species
        outputs.build = inputs.tcount.output.build
        outputs.source = inputs.source
        outputs.feature_type = "gene"

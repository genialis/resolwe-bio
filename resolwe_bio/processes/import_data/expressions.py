"""Upload Expressions."""

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
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


def parse_expression_file(exp_file, exp_type):
    """Parse expression file to a Pandas dataframe."""
    with gzip.open(exp_file) as exp:
        df = pd.read_csv(exp, sep="\t", float_precision="round_trip")

        df.rename(
            index=str,
            columns={
                "Gene": "FEATURE_ID",
                "Expression": exp_type,
            },
            inplace=True,
        )
        # Cast FEATURE_ID column to string
        df["FEATURE_ID"] = df["FEATURE_ID"].astype("str")
        # Remove any possible empty rows from the input file
        df.dropna(inplace=True)

    return df


def prepare_expression_set(exp, exp_type, feature_dict, outfile_name, rc=None):
    """Prepare expression set output data."""
    exp = parse_expression_file(exp_file=exp, exp_type=exp_type)
    exp["GENE_SYMBOL"] = exp["FEATURE_ID"].map(feature_dict)
    input_features = exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    # Merge expression values and reorder columns
    if rc:
        rc = parse_expression_file(exp_file=rc, exp_type="RAW_COUNT")
        exp_set = exp.merge(rc, on="FEATURE_ID")
        columns = ["FEATURE_ID", "GENE_SYMBOL", "RAW_COUNT", exp_type]
    else:
        exp_set = exp
        columns = ["FEATURE_ID", "GENE_SYMBOL", exp_type]
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


def expression_to_storage(infile, outfile):
    """Convert expressions file to JSON format."""

    def isfloat(value):
        """Check if value is float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    with io.TextIOWrapper(io.BufferedReader(gzip.open(infile))) as f:
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

    with open(file=outfile, mode="wt") as f:
        json.dump(exp, f)

    return outfile


def replace_extension(infile):
    """Replace extensions of file."""
    extensions = "".join(Path(str(infile)).suffixes[-2:])
    new_ext = ".tab.gz"
    outfile = str(infile).replace(extensions, new_ext)
    return outfile


class UploadExpression(ProcessBio):
    """Upload expression data.

    Upload expression data by providing raw expression data (read counts)
    and/or normalized expression data together with the associated data
    normalization type.
    """

    slug = "upload-expression"
    name = "Expression data"
    process_type = "data:expression"
    version = "2.6.0"
    category = "Import"
    data_name = "{{ exp_name }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "network": True,
        },
    }

    class Input:
        """Input fields to process UploadExpression."""

        rc = FileField(
            label="Read counts (raw expression)",
            description="Reads mapped to genomic features (raw count data). "
            "Supported extensions: .txt.gz (preferred), .tab.*, .txt.* or .tsv.*",
            required=False,
        )
        exp = FileField(
            label="Normalized expression",
            description="Normalized expression data. Supported extensions: .tab.gz "
            "(preferred), .tab.*, .txt.* or .tsv.*",
            required=False,
        )
        exp_name = StringField(
            label="Expression name",
        )
        exp_type = StringField(
            label="Normalization type",
            description="Normalization type",
            required=False,
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
        build = StringField(
            label="Build", description="Genome build or annotation version."
        )
        feature_type = StringField(
            label="Feature type",
            allow_custom_choice=True,
            default="gene",
            choices=[
                ("gene", "gene"),
                ("transcript", "transcript"),
                ("exon", "exon"),
            ],
        )

    class Output:
        """Output fields to process UploadExpression."""

        exp = FileField(label="Normalized expression")
        rc = FileField(
            label="Read counts",
            required=False,
            description="Reads mapped to genomic features.",
        )
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run analysis."""

        supported_extensions = (".txt", ".tab", ".tsv")

        if not inputs.exp and not inputs.rc:
            self.error("Please provide raw or/and normalized expression files.")

        elif inputs.exp and not inputs.exp_type:
            self.error(
                "Please provide normalization type together with normalized expressions."
            )

        elif not inputs.exp and inputs.exp_type and inputs.rc:
            self.error("Please provide raw or/and normalized expression files.")

        elif inputs.rc and not inputs.exp and not inputs.exp_type:
            rc = inputs.rc.import_file(imported_format="compressed")
            exp = inputs.rc.import_file(imported_format="compressed")
            exp_type = "RAW_COUNT"
            stem = Path(rc).stem

        elif inputs.exp and inputs.exp_type and not inputs.rc:
            exp = inputs.exp.import_file(imported_format="compressed")
            stem = Path(exp).stem
            exp_type = inputs.exp_type

        else:
            rc = inputs.rc.import_file(imported_format="compressed")
            exp = inputs.exp.import_file(imported_format="compressed")
            stem = Path(rc).stem
            stem_exp = Path(exp).stem
            if not stem_exp.endswith(supported_extensions):
                self.error(
                    f"The imported file has unsupported file name extension. "
                    f"The supported extensions are {supported_extensions}."
                )
            exp_type = inputs.exp_type

        if not stem.endswith(supported_extensions):
            self.error(
                "The imported file has unsupported file name extension. "
                f"The supported extensions are {supported_extensions}."
            )
        name = stem[:-4]

        # Save the abundance estimates to JSON storage
        expression_to_storage(infile=exp, outfile="json.txt")

        # Prepare the expression set outputs
        feature_ids = pd.read_csv(exp, sep="\t", index_col="Gene").index.tolist()

        feature_filters = {
            "source": inputs.source,
            "species": inputs.species,
            "feature_id__in": feature_ids,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        if inputs.rc and inputs.exp:
            prepare_expression_set(
                exp=exp,
                exp_type=exp_type,
                feature_dict=feature_ids_to_names,
                outfile_name=f"{name}_expressions",
                rc=rc,
            )
        else:
            prepare_expression_set(
                exp=exp,
                exp_type=exp_type,
                feature_dict=feature_ids_to_names,
                outfile_name=f"{name}_expressions",
            )

        # Change suffixes of exp file
        exp_final = replace_extension(infile=exp)
        Path(exp).rename(exp_final)
        exp = Path(exp_final).name

        if inputs.rc and inputs.exp:
            # Change suffixes of rc file
            rc_final = replace_extension(infile=rc)
            Path(rc).rename(rc_final)
            rc = Path(rc_final).name
            outputs.rc = rc
        elif inputs.rc and not inputs.exp:
            rc = exp
            outputs.rc = rc

        outputs.exp_type = exp_type
        outputs.exp = exp
        outputs.exp_json = "json.txt"
        outputs.exp_set = f"{name}_expressions.txt.gz"
        outputs.exp_set_json = f"{name}_expressions.json"
        outputs.source = inputs.source
        outputs.species = inputs.species
        outputs.build = inputs.build
        outputs.feature_type = inputs.feature_type


class UploadExpressionCuffnorm(ProcessBio):
    """Upload expression data by providing Cuffnorm results."""

    slug = "upload-expression-cuffnorm"
    name = "Expression data (Cuffnorm)"
    process_type = "data:expression"
    version = "1.8.0"
    category = "Import"
    data_name = '{{ exp.file|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "network": True,
        },
    }

    class Input:
        """Input fields for UploadExpressionCuffnorm."""

        exp = FileField(label="Normalized expression")
        cxb = DataField(
            "cufflinks:cuffquant",
            label="Cuffquant analysis",
            description="Cuffquant analysis.",
        )
        exp_type = StringField(
            label="Normalization type",
            default="Cuffnorm",
        )

    class Output:
        """Output fields for UploadExpressionCuffnorm."""

        exp = FileField(
            label="Normalized expression",
            description="Normalized expression",
        )
        exp_json = JsonField(
            label="Expression (json)",
        )
        exp_type = StringField(
            label="Expression type",
        )
        exp_set = FileField(
            label="Expressions",
        )
        exp_set_json = JsonField(
            label="Expressions (json)",
        )
        source = StringField(
            label="Gene ID source",
        )
        species = StringField(label="Species")
        build = StringField(
            label="Build",
        )
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run analysis."""

        if inputs.exp and not inputs.exp_type:
            self.error(
                "Please provide normalization type together with normalized expressions."
            )

        elif inputs.exp and inputs.exp_type and inputs.cxb:
            exp = inputs.exp.import_file(imported_format="compressed")
            stem = Path(exp).stem
            name = stem[:-4]

            # Save the abundance estimates to JSON storage
            expression_to_storage(infile=exp, outfile="json.txt")

            # Prepare the expression set outputs
            feature_ids = pd.read_csv(exp, sep="\t", index_col="Gene").index.tolist()

            feature_filters = {
                "source": inputs.cxb.output.source,
                "species": inputs.cxb.output.species,
                "feature_id__in": feature_ids,
            }

            feature_ids_to_names = {
                f.feature_id: f.name for f in self.feature.filter(**feature_filters)
            }

            prepare_expression_set(
                exp=exp,
                exp_type=inputs.exp_type,
                feature_dict=feature_ids_to_names,
                outfile_name=f"{name}_expressions",
            )

        outputs.exp_type = inputs.exp_type
        outputs.exp = exp
        outputs.exp_json = "json.txt"
        outputs.exp_set = f"{name}_expressions.txt.gz"
        outputs.exp_set_json = f"{name}_expressions.json"
        outputs.source = inputs.cxb.output.source
        outputs.species = inputs.cxb.output.species
        outputs.build = inputs.cxb.output.build
        outputs.feature_type = "gene"

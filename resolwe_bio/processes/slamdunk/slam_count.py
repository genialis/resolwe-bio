"""Convert collapsed Slamdunk tcount data to the expressions data type."""
import os

import pandas as pd
from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, JsonField, Process, StringField


def prepare_expressions(infile, rc_file, tpm_file):
    """Compute TPM values from T>C read counts."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        usecols=["gene_name", "length", "tcReadCount"],
        index_col="gene_name",
        dtype={"gene_name": str, "length": int, "tcReadCount": int,},
    )
    exp["rpk"] = exp.apply(lambda x: (x.tcReadCount * 1e3 / x.length), axis=1)
    rpk_sum = exp["rpk"].sum()
    exp["tpm"] = exp.apply(lambda x: (x.rpk / rpk_sum * 1e6), axis=1)
    rc = exp["tcReadCount"].to_csv(
        rc_file, index_label="Gene", header=["Expression"], sep="\t", compression="gzip"
    )
    tpm = exp["tpm"].to_csv(
        tpm_file,
        index_label="Gene",
        header=["Expression"],
        sep="\t",
        compression="gzip",
    )
    return (rc, tpm)


class SlamCount(Process):
    """Convert gene-level Slamdunk data to the expression data type."""

    slug = "slam-count"
    process_type = "data:expression:slamdunk"
    name = "Slam count"
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/slamdunk:1.0.0"},},
        "resources": {"cores": 1, "memory": 4096, "network": True,},
    }
    category = "Slamdunk"
    entity = {
        "type": "sample",
    }
    data_name = '{{ tcount|sample_name|default("?") }}'
    version = "1.0.2"

    class Input:
        """Input fields for SlamCount."""

        tcount = DataField("alleyoop:collapse", label="Collapsed Slamdunk results")
        source = StringField(
            label="Gene ID source",
            default="ENSEMBL",
            choices=[("ENSEMBL", "ENSEMBL"), ("UCSC", "UCSC"),],
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
        basename = os.path.basename(inputs.tcount.tcount.path)
        assert basename.endswith(".txt")
        name = basename[:-4]

        rc_file = name + "_rc.txt.gz"
        tpm_file = name + "_tmp.txt.gz"

        prepare_expressions(inputs.tcount.tcount.path, rc_file, tpm_file)

        for exp_file in [rc_file, tpm_file]:
            if not os.path.isfile(exp_file):
                self.error(
                    "Failed to parse tcout file. {} file was not created".format(
                        exp_file
                    )
                )

        # Save the abundance estimates to JSON storage
        Cmd["expression2storage.py"]("--output", "json.txt", tpm_file)

        # Prepare expression set file with feature_id -> gene_id mappings
        exp_set_args = [
            "--expressions",
            rc_file,
            "--source_db",
            inputs.source,
            "--species",
            inputs.tcount.species,
            "--output_name",
            name + "_expressions",
            "--norm_expressions",
            tpm_file,
            "--norm_expressions_type",
            "TPM",
        ]
        return_code, _, _ = Cmd["create_expression_set.py"][exp_set_args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("Error while preparing the expression set file.")

        outputs.exp = tpm_file
        outputs.exp_json = "json.txt"
        outputs.exp_type = "TPM"
        outputs.rc = rc_file
        outputs.exp_set = name + "_expressions.txt.gz"
        outputs.exp_set_json = name + "_expressions.json"
        outputs.species = inputs.tcount.species
        outputs.build = inputs.tcount.build
        outputs.source = inputs.source
        outputs.feature_type = "gene"

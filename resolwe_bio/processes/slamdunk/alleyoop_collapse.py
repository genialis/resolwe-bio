"""Run Alleyoop collapse tool on Slamdunk results."""
import os

import pandas as pd
from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, StringField

from resolwe_bio.process.runtime import ProcessBio


def compute_tpm(tcount):
    """Normalize readCount column to TPM values."""
    exp = pd.read_csv(
        tcount,
        sep="\t",
        index_col="gene_name",
    )
    exp["rpk"] = exp.apply(lambda x: (x.readCount * 1e3 / x.length), axis=1)
    rpk_sum = exp["rpk"].sum()
    exp["readsTPM"] = exp.apply(lambda x: (x.rpk / rpk_sum * 1e6), axis=1)
    return exp


class AlleyoopCollapse(ProcessBio):
    """Run Alleyoop collapse tool on Slamdunk results."""

    slug = "alleyoop-collapse"
    process_type = "data:alleyoop:collapse"
    name = "Alleyoop collapse"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/slamdunk:2.0.0"},
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
            "network": True,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "Slamdunk"
    data_name = '{{ slamdunk|sample_name|default("?") }}'
    version = "1.3.0"

    class Input:
        """Input fields for SlamdunkAllPaired."""

        slamdunk = DataField("alignment:bam:slamdunk", label="Slamdunk results")
        source = StringField(
            label="Gene ID source",
            default="ENSEMBL",
            choices=[
                ("ENSEMBL", "ENSEMBL"),
                ("UCSC", "UCSC"),
            ],
        )

    class Output:
        """Output fields to process SlamdunkAllPaired."""

        tcount = FileField(label="Count report containing SLAMSeq statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.slamdunk.output.tcount.path)
        assert basename.endswith(".tsv")
        name = basename[:-4]

        args = [
            "-o",
            ".",
            "-t",
            self.requirements.resources.cores,
        ]

        return_code, _, _ = Cmd["alleyoop"]["collapse"][args][
            inputs.slamdunk.output.tcount.path
        ] & TEE(retcode=None)
        if return_code:
            self.error("Alleyoop collapse analysis failed.")

        collapsed_output = name + "_collapsed.txt"
        os.rename(name + "_collapsed.csv", collapsed_output)

        # normalize to TPM
        tcount_tpm = compute_tpm(collapsed_output)

        # Map gene symbols to feature IDs
        feature_dict = {}
        out_columns = [
            "gene_symbol",
            "length",
            "readsCPM",
            "readsTPM",
            "conversionRate",
            "Tcontent",
            "coverageOnTs",
            "conversionsOnTs",
            "readCount",
            "tcReadCount",
            "multimapCount",
        ]

        input_features = tcount_tpm.index.tolist()

        features = {
            "source": inputs.source,
            "species": inputs.slamdunk.output.species,
            "feature_id__in": input_features,
        }
        feature_dict = {f.feature_id: f.name for f in self.feature.filter(**features)}

        tcount_tpm["gene_symbol"] = tcount_tpm.index.map(feature_dict)
        tcount_tpm.to_csv(collapsed_output, columns=out_columns, sep="\t")

        outputs.tcount = collapsed_output
        outputs.species = inputs.slamdunk.output.species
        outputs.build = inputs.slamdunk.output.build

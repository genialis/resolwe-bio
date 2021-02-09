"""Run Alleyoop summary tool on Slamdunk results."""
import os

import pandas as pd
from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, ListField, Process, StringField


def change_filename_column(infile, outfile):
    """Change FileName column in alleyoop summary output."""
    with open(infile, "r") as summary, open(outfile, "a") as out:
        header = summary.readline()
        report = pd.read_csv(summary, sep="\t", na_values=[], keep_default_na=False)
        report["FileName"] = report["FileName"].apply(os.path.basename)
        out.write(header)
        report.to_csv(out, sep="\t", index=False)


class AlleyoopSummary(Process):
    """Run Alleyoop summary."""

    slug = "alleyoop-summary"
    process_type = "data:alleyoop:summary"
    name = "Alleyoop summary"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/slamdunk:2.0.0"},
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "Slamdunk"
    data_name = '{{ slamdunk.0|sample_name|default("?") }}'
    version = "1.1.1"

    class Input:
        """Input fields for AlleyoopSummary."""

        slamdunk = ListField(
            DataField(
                data_type="alignment:bam:slamdunk",
                description="Select one or multiple data objects from slamdunk process.",
            ),
            label="Slamdunk results",
        )

    class Output:
        """Output fields to process AlleyoopSummary."""

        report = FileField(label="Tab-separated file with mapping statistics")
        plot_data = FileField(
            label="PCA values of the samples based on T>C read counts in regions of interest.",
            required=False,
        )
        plot = FileField(
            label="PCA plot of the samples based on T>C read counts in regions of interest.",
            required=False,
        )
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        for dunk in inputs.slamdunk:
            basename = os.path.basename(dunk.output.bam.path)
            assert basename.endswith(".bam")
            name = basename[:-4]
            (Cmd["ln"]["-s", dunk.output.tcount.path, f"{name}_tcount.tsv"])()

        report_file = "summary.txt"
        args = [
            "-o",
            report_file,
            "-t",
            ".",
        ]

        bam_paths = [dunk.output.bam.path for dunk in inputs.slamdunk]

        return_code, _, _ = Cmd["alleyoop"]["summary"][args][bam_paths] & TEE(
            retcode=None
        )
        if return_code:
            self.error("Alleyoop summary analysis failed.")

        if len(inputs.slamdunk) > 1:
            pca_data = "alleyoop_summary_PCA.txt"
            pca_plot = "alleyoop_summary_PCA.pdf"
            if os.path.isfile(pca_data) and os.path.isfile(pca_plot):
                outputs.plot_data = pca_data
                outputs.plot = pca_plot
            else:
                self.error("Failed to create PCA plot.")

        out_file = "alleyoop_summary.txt"
        change_filename_column(report_file, out_file)

        outputs.report = out_file
        outputs.species = inputs.slamdunk[0].output.species
        outputs.build = inputs.slamdunk[0].output.build

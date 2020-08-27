"""Run Alleyoop rates tool on Slamdunk results."""
import os

from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, Process, StringField


class AlleyoopRates(Process):
    """Run Alleyoop rates."""

    slug = "alleyoop-rates"
    process_type = "data:alleyoop:rates"
    name = "Alleyoop rates"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "resolwebio/slamdunk:1.0.0"},
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
    data_name = '{{ slamdunk|sample_name|default("?") }}'
    version = "1.0.0"

    class Input:
        """Input fields for AlleyoopRates."""

        ref_seq = DataField(
            "seq:nucleotide", label="FASTA file containig sequences for aligning"
        )
        slamdunk = DataField("alignment:bam:slamdunk", label="Slamdunk results")

    class Output:
        """Output fields to process AlleyoopRates."""

        report = FileField(
            label="Tab-separated file containing the overall conversion rates"
        )
        plot = FileField(label="Overall conversion rate plot file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.slamdunk.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        args = [
            "-o",
            "rates",
            "-r",
            inputs.ref_seq.fasta.path,
        ]

        return_code, _, _ = Cmd["alleyoop"]["rates"][args][
            inputs.slamdunk.bam.path
        ] & TEE(retcode=None)
        if return_code:
            self.error("Alleyoop rates analysis failed.")

        rates_file = os.path.join("rates", f"{name}_overallrates.csv")
        rates_file_renamed = os.path.join("rates", f"{name}_overallrates.txt")
        os.rename(rates_file, rates_file_renamed)

        outputs.report = rates_file_renamed
        outputs.plot = os.path.join("rates", f"{name}_overallrates.pdf")
        outputs.species = inputs.slamdunk.species
        outputs.build = inputs.slamdunk.build

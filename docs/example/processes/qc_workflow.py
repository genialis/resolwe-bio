"""Example QC workflow."""

from resolwe.process import (
    Data,
    DataField,
    Process,
)
from resolwe.process.models import Process as BioProcess


class ExampleQcWorkflow(Process):
    """Example workflow.

    This is an example workflow that takes FASTQ reads as input and performs
    quality control using FastQC tool.
    """

    slug = "workflow-qc-docs"
    name = "Example QC workflow"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    entity = {
        "type": "sample",
    }
    version = "1.0.0"
    process_type = "data:workflow:qc"
    category = "Pipeline"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Select sample(s) (FASTQ)",
            description="Reads in FASTQ file, single or paired end.",
        )

    class Output:
        """Output fields."""

        # workflows do not have output fields.

    def run(self, inputs, outputs):
        """Run the workflow."""

        # prepares inputs for the fastqc-paired-end process
        input_fastqc = {
            "reads": inputs.reads,
        }

        # trigger the FastQC analysis
        Data.create(
            process=BioProcess.get_latest("fastqc-paired-end"),
            input=input_fastqc,
            name=f"FastQC report ({inputs.reads.name})",
        )

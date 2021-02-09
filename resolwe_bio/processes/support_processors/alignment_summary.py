"""Produce a summary of alignment metrics."""
import gzip
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class AlignmentSummary(Process):
    """Produce a summary of alignment metrics from BAM file.

    Tool from Picard, wrapped by GATK4. See GATK
    CollectAlignmentSummaryMetrics for more information.
    """

    slug = "alignment-summary"
    name = "Picard AlignmentSummary"
    category = "Picard"
    process_type = "data:picard:summary"
    version = "2.1.1"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:5.2.0"}
        },
    }
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields for AlignmentSummary."""

        bam = DataField("alignment:bam", label="Alignment BAM file")
        genome = DataField("seq:nucleotide", label="Genome")

        adapters = DataField(
            "seq:nucleotide", label="Adapter sequences", required=False
        )

        validation_stringency = StringField(
            label="Validation stringency",
            description="Validation stringency for all SAM files read by this "
            "program. Setting stringency to SILENT can improve "
            "performance when processing a BAM file in which "
            "variable-length data (read, qualities, tags) do not "
            "otherwise need to be decoded. Default is STRICT.",
            choices=[
                ("STRICT", "STRICT"),
                ("LENIENT", "LENIENT"),
                ("SILENT", "SILENT"),
            ],
            default="STRICT",
        )

        insert_size = IntegerField(label="Maximum insert size", default=100000)

        pair_orientation = StringField(
            label="Pair orientation",
            default="null",
            choices=[
                ("null", "Unspecified"),
                ("FR", "FR"),
                ("RF", "RF"),
                ("TANDEM", "TANDEM"),
            ],
        )

        bisulfite = BooleanField(
            label="BAM file consists of bisulfite sequenced reads", default=False
        )

        assume_sorted = BooleanField(
            label="Sorted BAM file",
            description="If true the sort order in the header file will be ignored.",
            default=False,
        )

    class Output:
        """Output fields for AlignmentSummary."""

        report = FileField(label="Alignement metrics")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def get_sequences(self, fname):
        """Get a list of sequences from FASTA file."""
        if fname.endswith(".gz"):
            with gzip.open(fname, "rt") as fasta:
                return [line.strip() for line in fasta if not line.startswith(">")]
        else:
            with open(fname, "r") as fasta:
                return [line.strip() for line in fasta if not line.startswith(">")]

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.bam.output.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        metrics_file = f"{name}_alignment_metrics.txt"

        args = [
            "--INPUT",
            inputs.bam.output.bam.path,
            "--OUTPUT",
            metrics_file,
            "--REFERENCE_SEQUENCE",
            inputs.genome.output.fasta.path,
            "--VALIDATION_STRINGENCY",
            inputs.validation_stringency,
            "--MAX_INSERT_SIZE",
            inputs.insert_size,
            "--ASSUME_SORTED",
            inputs.assume_sorted,
            "--EXPECTED_PAIR_ORIENTATIONS",
            inputs.pair_orientation,
            "--IS_BISULFITE_SEQUENCED",
            inputs.bisulfite,
        ]

        if inputs.adapters:
            adapters_list = self.get_sequences(inputs.adapters.output.fasta.path)
            args.extend(["--ADAPTER_SEQUENCE", [adapters_list]])
        else:
            # Clear the default adapter list implemented in Picard.
            args.extend(["--ADAPTER_SEQUENCE", "null"])

        return_code, _, _ = Cmd["gatk"]["CollectAlignmentSummaryMetrics"][args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("CollectAlignmentSummaryMetrics tool failed.")

        outputs.report = metrics_file
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build

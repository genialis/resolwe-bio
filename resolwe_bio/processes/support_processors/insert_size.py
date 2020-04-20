"""Collect metrics about the insert size of a paired-end library."""
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    Process,
    SchedulingClass,
    StringField,
)


class InsertSizeMetrics(Process):
    """Collect metrics about the insert size of a paired-end library.

    Tool from Picard, wrapped by GATK4. See GATK
    CollectInsertSizeMetrics for more information.
    """

    slug = "insert-size"
    name = "Picard InsertSizeMetrics"
    category = "Picard"
    process_type = "data:picard:insert"
    version = "2.0.0"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/dnaseq:4.2.0"}},
    }
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields for InsertSizeMetrics."""

        bam = DataField("alignment:bam", label="Alignment BAM file")
        genome = DataField("seq:nucleotide", label="Genome")

        minimum_fraction = FloatField(
            label="Minimum fraction of reads in a category to be considered ",
            description="When generating the histogram, discard any data "
            "categories (out of FR, TANDEM, RF) that have fewer than this "
            "fraction of overall reads (Range: 0 and 0.5).",
            default=0.05,
        )

        include_duplicates = BooleanField(
            label="Include reads marked as duplicates in the insert size histogram",
            default=False,
        )

        deviations = FloatField(
            label="Deviations limit",
            description="Generate mean, standard deviation and plots by trimming "
            "the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. "
            "This is done because insert size data typically includes enough "
            "anomalous values from chimeras and other artifacts to make the "
            "mean and standard deviation grossly misleading regarding the real "
            "distribution.",
            default=10.0,
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

        assume_sorted = BooleanField(
            label="Sorted BAM file",
            description="If True, the sort order in the header file will be ignored.",
            default=False,
        )

    class Output:
        """Output fields for InsertSizeMetrics."""

        report = FileField(label="Insert size metrics")
        plot = FileField(label="Insert size histogram")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.bam.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        metrics_file = f"{name}_insert_size_metrics.txt"
        histogram_file = f"{name}_insert_size.pdf"

        args = [
            "--INPUT",
            inputs.bam.bam.path,
            "--OUTPUT",
            metrics_file,
            "--Histogram_FILE",
            histogram_file,
            "--REFERENCE_SEQUENCE",
            inputs.genome.fasta.path,
            "--DEVIATIONS",
            inputs.deviations,
            "--INCLUDE_DUPLICATES",
            inputs.include_duplicates,
            "--VALIDATION_STRINGENCY",
            inputs.validation_stringency,
            "--ASSUME_SORTED",
            inputs.assume_sorted,
        ]

        if 0 <= inputs.minimum_fraction <= 0.5:
            args.extend(["--MINIMUM_PCT", inputs.minimum_fraction])
        else:
            self.warning(
                "Minimum fraction of reads should be between 0 and 0.5. "
                "Setting minimum fraction of reads to 0."
            )
            args.extend(["--MINIMUM_PCT", 0])

        return_code, _, _ = Cmd["gatk"]["CollectInsertSizeMetrics"][args] & TEE(
            retcode=None
        )
        if return_code:
            self.error("CollectInsertSizeMetrics tool failed.")

        outputs.report = metrics_file
        outputs.plot = histogram_file
        outputs.species = inputs.bam.species
        outputs.build = inputs.bam.build

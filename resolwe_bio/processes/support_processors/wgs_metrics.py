"""Collect metrics about coverage of whole genome sequencing."""
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


def replace_metrics_class(fname):
    """Replace metrics class name.

    This temporary fix is needed due to compatibility issue with GATK
    4.1.2.0 and MultiQC 1.8. MultiQC searches for CollectWgsMetrics
    instead of WgsMetrics in the report file. Note that this should be
    resolved in the 1.9 release of MultiQC.
    """
    with open(fname, "r") as report:
        newlines = []
        for line in report.readlines():
            if line == "## METRICS CLASS\tpicard.analysis.WgsMetrics\n":
                line = "## METRICS CLASS\tCollectWgsMetrics$WgsMetrics\n"
                newlines.append(line)
            else:
                newlines.append(line)
    with open(fname, "w") as report:
        for line in newlines:
            report.writelines(line)


class InsertSizeMetrics(Process):
    """Collect metrics about coverage of whole genome sequencing.

    Tool from Picard, wrapped by GATK4. See GATK
    CollectWgsMetrics for more information.
    """

    slug = "wgs-metrics"
    name = "Picard WGS Metrics"
    category = "Picard"
    process_type = "data:picard:wgsmetrics"
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
        """Input fields for CollectWgsMetrics."""

        bam = DataField("alignment:bam", label="Alignment BAM file")
        genome = DataField("seq:nucleotide", label="Genome")

        read_length = IntegerField(label="Average read length", default=150)

        create_histogram = BooleanField(
            label="Include data for base quality histogram in the metrics file",
            default=False,
        )

        advanced = BooleanField(
            label="Show advanced options",
            description="Inspect and modify parameters.",
            default=False,
        )

        class Options:
            """Options."""

            min_map_quality = IntegerField(
                label="Minimum mapping quality for a read to contribute coverage",
                default=20,
            )

            min_quality = IntegerField(
                label="Minimum base quality for a base to contribute coverage",
                description="N bases will be treated as having a base quality of "
                "negative infinity and will therefore be excluded from coverage "
                "regardless of the value of this parameter.",
                default=20,
            )

            coverage_cap = IntegerField(
                label="Maximum coverage cap",
                description="Treat positions with coverage exceeding this value as "
                "if they had coverage at this set value.",
                default=250,
            )

            accumulation_cap = IntegerField(
                label="Ignore positions with coverage above this value",
                description="At positions with coverage exceeding this value, "
                "completely ignore reads that accumulate beyond this value",
                default=100000,
            )

            count_unpaired = BooleanField(
                label="Count unpaired reads and paired reads with one end unmapped",
                default=False,
            )

            sample_size = IntegerField(
                label="Sample Size used for Theoretical Het Sensitivity sampling",
                default=10000,
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

        options = GroupField(Options, label="Options", hidden="!advanced")

    class Output:
        """Output fields for CollectWgsMetrics."""

        report = FileField(label="WGS metrics report")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.bam.output.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        metrics_file = f"{name}_wgs_metrics.txt"

        args = [
            "--INPUT",
            inputs.bam.output.bam.path,
            "--OUTPUT",
            metrics_file,
            "--REFERENCE_SEQUENCE",
            inputs.genome.output.fasta.path,
            "--READ_LENGTH",
            inputs.read_length,
            "--INCLUDE_BQ_HISTOGRAM",
            inputs.create_histogram,
            "--MINIMUM_MAPPING_QUALITY",
            inputs.options.min_map_quality,
            "--MINIMUM_BASE_QUALITY",
            inputs.options.min_quality,
            "--COVERAGE_CAP",
            inputs.options.coverage_cap,
            "--LOCUS_ACCUMULATION_CAP",
            inputs.options.accumulation_cap,
            "--COUNT_UNPAIRED",
            inputs.options.count_unpaired,
            "--SAMPLE_SIZE",
            inputs.options.sample_size,
            "--VALIDATION_STRINGENCY",
            inputs.options.validation_stringency,
        ]

        return_code, _, _ = Cmd["gatk"]["CollectWgsMetrics"][args] & TEE(retcode=None)
        if return_code:
            self.error("CollectWgsMetrics tool failed.")

        replace_metrics_class(metrics_file)

        outputs.report = metrics_file
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build

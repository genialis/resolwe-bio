"""Produce metrics for RRBS data based on the methylation status ."""
import os

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class CollectRrbsMetrics(Process):
    """Produce metrics for RRBS data based on the methylation status.

    This tool uses reduced representation bisulfite sequencing (Rrbs)
    data to determine cytosine methylation status across all reads of
    a genomic DNA sequence.

    Tool is wrapped by GATK4. See GATK
    CollectRrbsMetrics for more information.
    """

    slug = "rrbs-metrics"
    name = "Picard CollectRrbsMetrics"
    category = "Picard"
    process_type = "data:picard:rrbs"
    version = "2.0.0"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/dnaseq:4.2.0"}},
        "resources": {
            "memory": 32768,
        },
    }
    data_name = '{{ bam|sample_name|default("?") }}'

    class Input:
        """Input fields for CollectRrbsMetrics."""

        bam = DataField("alignment:bam", label="Alignment BAM file")
        genome = DataField("seq:nucleotide", label="Genome")

        min_quality = IntegerField(
            label="Threshold for base quality of a C base before it is considered",
            default=20,
        )

        next_base_quality = IntegerField(
            label="Threshold for quality of a base next to a C before the C base is considered",
            default=10,
        )

        min_lenght = IntegerField(label="Minimum read length", default=5)

        mismatch_rate = FloatField(
            label="Maximum fraction of mismatches in a read to be considered (Range: 0 and 1)",
            default=0.1,
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
            description="If true the sort order in the header file will be ignored.",
            default=False,
        )

    class Output:
        """Output fields for CollectRrbsMetrics."""

        report = FileField(label="RRBS summary metrics")
        detailed_report = FileField(label="Detailed RRBS report")
        plot = FileField(label="QC plots")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.bam.bam.path)
        assert basename.endswith(".bam")
        name = basename[:-4]

        args = [
            "--INPUT",
            inputs.bam.bam.path,
            "--REFERENCE",
            inputs.genome.fasta.path,
            "--METRICS_FILE_PREFIX",
            name,
            "--C_QUALITY_THRESHOLD",
            inputs.min_quality,
            "--NEXT_BASE_QUALITY_THRESHOLD",
            inputs.next_base_quality,
            "--MINIMUM_READ_LENGTH",
            inputs.min_lenght,
            "--VALIDATION_STRINGENCY",
            inputs.validation_stringency,
            "--ASSUME_SORTED",
            inputs.assume_sorted,
        ]

        if 0 <= inputs.mismatch_rate <= 1:
            args.extend(["--MAX_MISMATCH_RATE", inputs.mismatch_rate])

        return_code, _, _ = Cmd["gatk"]["CollectRrbsMetrics"][args] & TEE(retcode=None)
        if return_code:
            self.error("CollectRrbsMetrics tool failed.")

        report_file = f"{name}_rrbs_summary_metrics.txt"
        os.rename(f"{name}.rrbs_summary_metrics", report_file)

        detailed_file = f"{name}_rrbs_detail_metrics.txt"
        os.rename(f"{name}.rrbs_detail_metrics", detailed_file)

        out_plot = f"{name}_rrbs_qc.pdf"
        os.rename(f"{name}.rrbs_qc.pdf", out_plot)

        outputs.report = report_file
        outputs.detailed_report = detailed_file
        outputs.plot = out_plot
        outputs.species = inputs.bam.species
        outputs.build = inputs.bam.build

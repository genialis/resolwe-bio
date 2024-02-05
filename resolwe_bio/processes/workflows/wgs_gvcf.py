"""WGS primary analysis (GVCF)."""

from resolwe.process import (
    BooleanField,
    Data,
    DataField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    StringField,
)
from resolwe.process.models import Process as BioProcess


class WorkflowWgsGvcf(Process):
    """Whole genome sequencing pipeline (GATK GVCF).

    The pipeline follows GATK best practices recommendations and prepares
    single-sample paired-end sequencing data for a joint-genotyping step.

    The pipeline steps include read trimming (Trimmomatic), read alignment
    (BWA-MEM2), marking of duplicates (Picard MarkDuplicates), recalibration
    of base quality scores (ApplyBQSR) and calling of variants
    (GATK HaplotypeCaller in GVCF mode). The QC reports (FASTQC report,
    Picard AlignmentSummaryMetrics, CollectWgsMetrics and InsertSizeMetrics)
    are summarized using MultiQC.
    """

    slug = "workflow-wgs-gvcf"
    name = "WGS analysis (GVCF)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1",
            },
        },
    }
    data_name = "WGS GVCF analysis ({{ reads|name|default('?') if reads else aligned_reads|name|default('?') }})"
    version = "2.3.0"
    process_type = "data:workflow:wgs:gvcf"
    category = "Pipeline"
    entity = {
        "type": "sample",
    }

    class Input:
        """Input fields."""

        reads = DataField(
            "reads:fastq:paired",
            label="Input sample (FASTQ)",
            required=False,
            disabled="aligned_reads",
            description="Input data in FASTQ format. This input type allows for optional "
            "read trimming procedure and is mutually exclusive with the BAM input file type.",
        )
        aligned_reads = DataField(
            "alignment:bam",
            label="Input sample (BAM)",
            required=False,
            disabled="reads",
            description="Input data in BAM format. This input file type is mutually exclusive "
            "with the FASTQ input file type and does not allow for read trimming procedure.",
        )
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")
        bwa_index = DataField("index:bwamem2", label="BWA genome index")
        known_sites = ListField(
            DataField("variants:vcf"), label="Known sites of variation (VCF)"
        )

        class Trimming:
            """Trimming parameters."""

            enable_trimming = BooleanField(
                label="Trim and quality filter input data",
                description="Enable or disable adapter trimming and QC filtering procedure.",
                default=False,
            )

            adapters = DataField(
                "seq:nucleotide",
                label="Adapter sequences",
                required=False,
                description="Adapter sequences in FASTA format that will "
                "be removed from the reads.",
                disabled="!trimming_options.enable_trimming",
            )

            seed_mismatches = IntegerField(
                label="Seed mismatches",
                required=False,
                disabled="!trimming_options.adapters",
                description="Specifies the maximum mismatch count which "
                "will still allow a full match to be performed. This field "
                "is required to perform adapter trimming.",
            )

            simple_clip_threshold = IntegerField(
                label="Simple clip threshold",
                required=False,
                disabled="!trimming_options.adapters",
                description="Specifies how accurate the match between any "
                "adapter sequence must be against a read. This field is "
                "required to perform adapter trimming.",
            )

            min_adapter_length = IntegerField(
                label="Minimum adapter length",
                default=8,
                disabled="!trimming_options.seed_mismatches && "
                "!trimming_options.simple_clip_threshold && "
                "!trimming_options.palindrome_clip_threshold",
                description="In addition to the alignment score, palindrome "
                "mode can verify that a minimum length of adapter has been "
                "detected. If unspecified, this defaults to 8 bases, for "
                "historical reasons. However, since palindrome mode has a "
                "very low false positive rate, this can be safely reduced, "
                "even down to 1, to allow shorter adapter fragments to be "
                "removed.",
            )

            palindrome_clip_threshold = IntegerField(
                label="Palindrome clip threshold",
                required=False,
                disabled="!trimming_options.adapters",
                description="Specifies how accurate the match between the "
                "two adapter ligated reads must be for PE palindrome read "
                "alignment. This field is required to perform adapter "
                "trimming.",
            )

            leading = IntegerField(
                label="Leading quality",
                required=False,
                description="Remove low quality bases from the beginning, "
                "if below a threshold quality.",
                disabled="!trimming_options.enable_trimming",
            )

            trailing = IntegerField(
                label="Trailing quality",
                required=False,
                description="Remove low quality bases from the end, if "
                "below a threshold quality.",
                disabled="!trimming_options.enable_trimming",
            )

            minlen = IntegerField(
                label="Minimum length",
                required=False,
                description="Drop the read if it is below a specified length.",
                disabled="!trimming_options.enable_trimming",
            )

        class GatkOptions:
            """Options."""

            intervals = DataField(
                "bed",
                label="Intervals BED file",
                description="Use intervals BED file to limit the analysis to "
                "the specified parts of the genome.",
                required=False,
            )

            contamination = IntegerField(
                label="Contamination fraction",
                default=0,
                description="Fraction of contamination in sequencing "
                "data (for all samples) to aggressively remove.",
            )

        class AlignmentSummary:
            """AlignmentSummary parameters."""

            adapters = DataField(
                "seq:nucleotide",
                label="Adapter sequences",
                required=False,
            )

            max_insert_size = IntegerField(
                label="Maximum insert size",
                default=100000,
            )

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

        class PicardWGSMetrics:
            """PicardWGSMetrics parameters."""

            read_length = IntegerField(
                label="Average read length",
                default=150,
            )

            min_map_quality = IntegerField(
                label="Minimum mapping quality for a read to contribute coverage",
                default=20,
            )

            min_quality = IntegerField(
                label="Minimum base quality for a base to contribute coverage",
                default=20,
                description="N bases will be treated as having a base quality of "
                "negative infinity and will therefore be excluded from "
                "coverage regardless of the value of this parameter.",
            )

            coverage_cap = IntegerField(
                label="Maximum coverage cap",
                default=250,
                description="Treat positions with coverage exceeding this "
                "value as if they had coverage at this set value.",
            )

            accumulation_cap = IntegerField(
                label="Ignore positions with coverage above this value",
                default=100000,
                description="At positions with coverage exceeding this value, "
                "completely ignore reads that accumulate beyond this value.",
            )

            sample_size = IntegerField(
                label="Sample size used for Theoretical Het Sensitivity sampling",
                default=10000,
            )

        class InsertSizeMetrics:
            """InsertSizeMetrics parameters."""

            minimum_fraction = FloatField(
                label="Minimum fraction of reads in a category to be considered",
                default=0.05,
                description="When generating the histogram, discard any data "
                "categories (out of FR, TANDEM, RF) that have fewer than "
                "this fraction of overall reads (Range: 0 and 0.5).",
            )

            include_duplicates = BooleanField(
                label="Include reads marked as duplicates in the insert size histogram",
                default=False,
            )

            deviations = FloatField(
                label="Deviations limit",
                default=10.0,
                description="Generate mean, standard deviation and plots "
                "by trimming the data down to MEDIAN + DEVIATIONS * "
                "MEDIAN_ABSOLUTE_DEVIATION. This is done because insert "
                "size data typically includes enough anomalous values "
                "from chimeras and other artifacts to make the mean and "
                "standard deviation grossly misleading regarding the real "
                "distribution.",
            )

        trimming_options = GroupField(Trimming, label="Trimming options")

        gatk_options = GroupField(GatkOptions, label="GATK options")

        alignment_summary = GroupField(
            AlignmentSummary, label="Alignment summary options"
        )

        wgs_metrics = GroupField(PicardWGSMetrics, label="Picard WGS metrics options")

        insert_size = GroupField(
            InsertSizeMetrics, label="Picard InsertSizeMetrics options"
        )

    class Output:
        """Output fields."""

    def run(self, inputs, outputs):
        """Run the workflow."""
        if not inputs.reads and not inputs.aligned_reads:
            self.error("Please provide FASTQ or BAM input files.")
        if inputs.reads and inputs.aligned_reads:
            self.error(
                "Please provide input data in either FASTQ or aligned BAM format, not both."
            )

        preprocess_inputs = {
            "ref_seq": inputs.ref_seq,
            "bwa_index": inputs.bwa_index,
            "known_sites": inputs.known_sites,
        }

        if inputs.reads and inputs.trimming_options.enable_trimming:
            trimmomatic = Data.create(
                process=BioProcess.get_latest(slug="trimmomatic-paired"),
                input={
                    "reads": inputs.reads,
                    "illuminaclip": {
                        "adapters": inputs.trimming_options.adapters,
                        "seed_mismatches": inputs.trimming_options.seed_mismatches,
                        "simple_clip_threshold": inputs.trimming_options.simple_clip_threshold,
                        "palindrome_clip_threshold": inputs.trimming_options.palindrome_clip_threshold,
                        "min_adapter_length": inputs.trimming_options.min_adapter_length,
                    },
                    "trim_bases": {
                        "trailing": inputs.trimming_options.trailing,
                        "leading": inputs.trimming_options.leading,
                    },
                    "reads_filtering": {"minlen": inputs.trimming_options.minlen},
                },
            )
            preprocess_inputs.update(reads=trimmomatic)
        elif inputs.reads and not inputs.trimming_options.enable_trimming:
            preprocess_inputs.update(reads=inputs.reads)
        else:
            preprocess_inputs.update(aligned_reads=inputs.aligned_reads)

        bam = Data.create(
            process=BioProcess.get_latest(slug="wgs-preprocess-bwa2"),
            input=preprocess_inputs,
        )

        Data.create(
            process=BioProcess.get_latest(slug="gatk-haplotypecaller-gvcf"),
            input={
                "bam": bam,
                "ref_seq": inputs.ref_seq,
                "options": {
                    "intervals": inputs.gatk_options.intervals,
                    "contamination": inputs.gatk_options.contamination,
                },
            },
        )

        alignment_summary_inputs = {
            "bam": bam,
            "genome": inputs.ref_seq,
            "insert_size": inputs.alignment_summary.max_insert_size,
            "pair_orientation": inputs.alignment_summary.pair_orientation,
            "bisulfite": False,
            "assume_sorted": True,
        }
        if inputs.alignment_summary.adapters:
            alignment_summary_inputs.update(
                {"adapters": inputs.alignment_summary.adapters}
            )

        summary = Data.create(
            process=BioProcess.get_latest(slug="alignment-summary"),
            input=alignment_summary_inputs,
        )

        wgs_metrics = Data.create(
            process=BioProcess.get_latest(slug="wgs-metrics"),
            input={
                "bam": bam,
                "genome": inputs.ref_seq,
                "read_length": inputs.wgs_metrics.read_length,
                "create_histogram": False,
                "options": {
                    "min_map_quality": inputs.wgs_metrics.min_map_quality,
                    "coverage_cap": inputs.wgs_metrics.coverage_cap,
                    "accumulation_cap": inputs.wgs_metrics.accumulation_cap,
                    "count_unpaired": False,
                    "sample_size": inputs.wgs_metrics.sample_size,
                },
            },
        )

        insert_size = Data.create(
            process=BioProcess.get_latest(slug="insert-size"),
            input={
                "bam": bam,
                "genome": inputs.ref_seq,
                "minimum_fraction": inputs.insert_size.minimum_fraction,
                "include_duplicates": inputs.insert_size.include_duplicates,
                "deviations": inputs.insert_size.deviations,
                "assume_sorted": True,
            },
        )
        multiqc_inputs = [
            bam,
            summary,
            wgs_metrics,
            insert_size,
        ]
        if inputs.reads:
            multiqc_inputs.append(inputs.reads)
        if inputs.reads and inputs.trimming_options.enable_trimming:
            multiqc_inputs.append(trimmomatic)
        Data.create(
            process=BioProcess.get_latest(slug="multiqc"),
            input={"data": multiqc_inputs},
        )

"""Preprocess data for WGS analysis using GATK best practices procedure."""
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


class WgsPreprocess(Process):
    """Prepare analysis ready BAM file.

    This process follows GATK best practices procedure to prepare
    analysis-ready BAM file. The steps included are read alignment using
    BWA MEM, marking of duplicates (Picard MarkDuplicates), BAM sorting,
    read-group assignment and base quality score recalibration (BQSR).
    """

    slug = "wgs-preprocess"
    name = "WGS preprocess data"
    process_type = "data:alignment:bam:wgs"
    version = "1.0.0"
    category = "GATK"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
            "storage": 400,
        },
    }
    data_name = '{{ reads|sample_name|default("?") }}'

    class Input:
        """Input fields to process WgsPreprocess."""

        reads = DataField("reads:fastq:paired", label="Input sample")
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")
        bwa_index = DataField("index:bwa", label="BWA genome index")
        known_sites = ListField(
            DataField("variants:vcf"), label="Known sites of variation (VCF)"
        )

        advanced = BooleanField(
            label="Show advanced options",
            description="Inspect and modify parameters.",
            default=False,
        )

        class AdvancedOptions:
            """Advanced options."""

            pixel_distance = IntegerField(
                label="--OPTICAL_DUPLICATE_PIXEL_DISTANCE",
                default=2500,
                description="Set the optical pixel distance, e.g. "
                "distance between clusters. Modify this parameter to "
                "ensure compatibility with older Illumina platforms.",
            )

        advanced_options = GroupField(
            AdvancedOptions, label="Advanced options", hidden="!advanced"
        )

    class Output:
        """Output fields to process WgsPreprocess."""

        bam = FileField(label="Analysis ready BAM file")
        bai = FileField(label="BAM file index")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")
        metrics_file = FileField(label="Metrics from MarkDuplicate process")

    def run(self, inputs, outputs):
        """Run analysis."""

        # Define output file names
        name = inputs.reads.entity_name
        if name == None or name == "":
            mate1_path = Path(inputs.reads.output.fastq[0].path).name
            assert mate1_path.endswith(".fastq.gz")
            name = mate1_path[:-9]

        index_fasta_name = Path(inputs.bwa_index.output.fasta.path).name
        bam = f"{name}.bam"
        bam_stats = f"{name}_stats.txt"
        metrics_file = f"{name}_markduplicates_metrics.txt"

        aligned_bam = "aligned.bam"
        marked_dups = "marked_duplicates.bam"
        sorted_bam = "sorted.bam"
        sorted_rg = "sorted_rg.bam"
        recal_table = "recal_data.csv"

        # Concatenate multi-lane read files
        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]]
            > "input_reads_mate1.fastq.gz"
        )()
        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq2]]
            > "input_reads_mate2.fastq.gz"
        )()

        self.progress(0.05)

        # Align reads with BWA MEM
        bwa_inputs = [
            "-K 100000000",
            "-v 3",
            f"-t {self.requirements.resources.cores}",
            "-Y",
            f"{Path(inputs.bwa_index.output.index.path) / index_fasta_name}",
            "input_reads_mate1.fastq.gz",
            "input_reads_mate2.fastq.gz",
        ]
        (
            Cmd["bwa"]["mem"][bwa_inputs] | Cmd["samtools"]["view"]["-1"]["-"]
            > aligned_bam
        )()

        self.progress(0.2)

        # Mark duplicates
        mark_duplicates_inputs = [
            "--INPUT",
            aligned_bam,
            "--OUTPUT",
            marked_dups,
            "--METRICS_FILE",
            metrics_file,
            "--VALIDATION_STRINGENCY",
            "SILENT",
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE",
            inputs.advanced_options.pixel_distance,
            "--ASSUME_SORT_ORDER",
            "queryname",
        ]

        return_code, _, _ = Cmd["gatk"]["MarkDuplicates"][mark_duplicates_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("MarkDuplicates analysis failed.")

        self.progress(0.3)

        # Sort BAM file and fix NM and UQ tags
        sort_inputs = [
            "--INPUT",
            marked_dups,
            "--OUTPUT",
            "/dev/stdout",
            "--SORT_ORDER",
            "coordinate",
            "--CREATE_INDEX",
            "false",
        ]
        set_tag_inputs = [
            "--INPUT",
            "/dev/stdin",
            "--OUTPUT",
            sorted_bam,
            "--CREATE_INDEX",
            "true",
            "--REFERENCE_SEQUENCE",
            inputs.ref_seq.output.fasta.path,
        ]
        (
            Cmd["gatk"]["SortSam"][sort_inputs]
            | Cmd["gatk"]["SetNmMdAndUqTags"][set_tag_inputs]
        )()

        self.progress(0.4)

        # Set the read group information (required by BaseRecalibrator)
        rg_inputs = [
            "--INPUT",
            sorted_bam,
            "--VALIDATION_STRINGENCY",
            "STRICT",
            "--OUTPUT",
            sorted_rg,
            "--RGLB",
            "WGS",
            "--RGPL",
            "ILLUMINA",
            "--RGPU",
            "X",
            "--RGSM",
            name,
        ]

        return_code, _, _ = Cmd["gatk"]["AddOrReplaceReadGroups"][rg_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("AddOrReplaceReadGroups tool failed.")

        self.progress(0.45)

        # BaseRecalibrator
        base_recal_inputs = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-I",
            sorted_rg,
            "--use-original-qualities",
            "-O",
            recal_table,
        ]
        # Add known sites to the input parameters of BaseRecalibrator.
        for site in inputs.known_sites:
            base_recal_inputs.extend(["--known-sites", f"{site.output.vcf.path}"])

        return_code, _, _ = Cmd["gatk"]["BaseRecalibrator"][base_recal_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("BaseRecalibrator analysis failed.")

        self.progress(0.6)

        # ApplyBQSR
        apply_bqsr_inputs = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-I",
            sorted_rg,
            "-O",
            bam,
            "-bqsr",
            recal_table,
            "--static-quantized-quals",
            "10",
            "--static-quantized-quals",
            "20",
            "--static-quantized-quals",
            "30",
            "--add-output-sam-program-record",
            "--use-original-qualities",
        ]
        return_code, _, _ = Cmd["gatk"]["ApplyBQSR"][apply_bqsr_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("ApplyBQSR analysis failed.")

        self.progress(0.9)

        # Index the BQSR BAM file
        return_code, _, _ = Cmd["samtools"]["index"][bam] & TEE(retcode=None)
        if return_code:
            self.error("Samtools index command failed.")

        (Cmd["samtools"]["flagstat"][f"{bam}"] > bam_stats)()
        self.progress(0.95)

        outputs.bam = bam
        outputs.bai = bam + ".bai"
        outputs.stats = bam_stats
        outputs.species = inputs.bwa_index.output.species
        outputs.build = inputs.bwa_index.output.build
        outputs.metrics_file = metrics_file

"""Preprocess data for WGS analysis using GATK best practices procedure."""

import os
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed, wrap_non_picklable_objects
from plumbum import TEE

from resolwe.process import (
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


def prepare_chromosome_sizes(fai_path, bed_path):
    """Prepare a BED file with chromosome sizes."""
    fai = pd.read_csv(
        fai_path,
        sep="\t",
        header=None,
        names=["chr", "length", "offset", "line_bases", "line_width"],
    )
    fai = fai.drop(columns=["offset", "line_bases", "line_width"])
    fai.insert(loc=1, column="start", value=0)
    fai.to_csv(bed_path, sep="\t", header=False, index=False)


def prepare_scattered_inputs(results_dir, pattern="*"):
    """Prepare the input arguments for scattered input files.

    This expects the files in results_dir to be named using four number
    interval notation used by GATK SplitIntervals (e.g. 0001-scattered).
    Names are used for sorting, which ensures the correct concatenation
    order.
    """
    input_list = []
    for scattered_output in sorted(results_dir.glob(pattern)):
        input_list.extend(["-I", scattered_output])
    return input_list


@delayed
@wrap_non_picklable_objects
def run_base_recalibration(
    intput_bam, known_sites, interval_path, ref_seq_path, tmp, parent_dir
):
    """Run BaseRecalibrator on a specifed interval."""
    recal_interval = f"{parent_dir.name}/{interval_path.stem}.recal_data.csv"

    base_recal_inputs = [
        "-R",
        ref_seq_path,
        "-I",
        intput_bam,
        "--use-original-qualities",
        "--tmp-dir",
        tmp,
        "-L",
        interval_path,
        "-O",
        recal_interval,
    ]
    # Add known sites to the input parameters of BaseRecalibrator.
    for site in known_sites:
        base_recal_inputs.extend(["--known-sites", f"{site.output.vcf.path}"])

    return_code, stdout, stderr = Cmd["gatk"]["BaseRecalibrator"][
        base_recal_inputs
    ] & TEE(retcode=None)
    if return_code:
        print(f"Error in {interval_path.stem} interval.", stdout, stderr)
    return return_code


@delayed
@wrap_non_picklable_objects
def run_apply_bqsr(
    intput_bam, recal_table, interval_path, ref_seq_path, tmp, parent_dir
):
    """Run ApplyBQSR on a specifed interval."""
    bqsr_interval_bam = f"{parent_dir.name}/{interval_path.stem}.bam"

    apply_bqsr_inputs = [
        "-R",
        ref_seq_path,
        "-I",
        intput_bam,
        "-O",
        bqsr_interval_bam,
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
        "-L",
        interval_path,
        "--tmp-dir",
        tmp,
    ]
    return_code, stdout, stderr = Cmd["gatk"]["ApplyBQSR"][apply_bqsr_inputs] & TEE(
        retcode=None
    )
    if return_code:
        print(f"Error in {interval_path.stem} interval.", stdout, stderr)
    return return_code


class WgsPreprocess_BWA2(Process):
    """Prepare analysis ready BAM file.

    This process follows GATK best practices procedure to prepare
    analysis-ready BAM file. The steps included are read alignment using
    BWA MEM2, marking of duplicates (Picard MarkDuplicates), BAM sorting,
    read-group assignment and base quality score recalibration (BQSR).
    """

    slug = "wgs-preprocess-bwa2"
    name = "WGS preprocess data with bwa-mem2"
    process_type = "data:alignment:bam:wgsbwa2"
    version = "1.4.0"
    category = "WGS"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
            "storage": 600,
        },
    }
    data_name = (
        "{{ reads|name|default('?') if reads else aligned_reads|name|default('?') }}"
    )

    class Input:
        """Input fields to process WgsPreprocess_BWA2."""

        reads = DataField(
            "reads:fastq:paired", label="Input sample (FASTQ)", required=False
        )
        aligned_reads = DataField(
            "alignment:bam", label="Input sample (BAM)", required=False
        )
        ref_seq = DataField("seq:nucleotide", label="Reference sequence")
        bwa_index = DataField("index:bwamem2", label="BWA-MEM2 genome index")
        known_sites = ListField(
            DataField("variants:vcf"), label="Known sites of variation (VCF)"
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

            n_jobs = IntegerField(
                label="Number of concurent jobs",
                description="Use a fixed number of jobs for quality score "
                "recalibration of determining it based on the number of "
                "available cores.",
                required=False,
            )

        advanced_options = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields to process WgsPreprocess_BWA2."""

        bam = FileField(label="Analysis ready BAM file")
        bai = FileField(label="BAM file index")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")
        metrics_file = FileField(label="Metrics from MarkDuplicate process")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        aligned_sam = "aligned.sam"
        aligned_bam = "aligned.bam"
        marked_dups = "marked_duplicates.bam"
        sorted_temp = "sorted_temp.bam"
        sorted_bam = "sorted.bam"
        sorted_rg = "sorted_rg.bam"
        recal_table = "recal_data.csv"

        index_fasta_name = Path(inputs.bwa_index.output.fasta.path).name

        if not inputs.reads and not inputs.aligned_reads:
            self.error("Please provide FASTQ or BAM input files.")
        if inputs.reads and inputs.aligned_reads:
            self.error(
                "Please provide input data in either FASTQ or aligned BAM format, not both."
            )

        if inputs.reads:
            # Define output file names
            name = inputs.reads.entity_name
            if not name:
                mate1_path = Path(inputs.reads.output.fastq[0].path).name
                assert mate1_path.endswith(".fastq.gz")
                name = mate1_path[:-9]

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

            # Align reads with BWA MEM2
            bwa_inputs = [
                "-K 100000000",
                "-v 3",
                f"-t {self.requirements.resources.cores}",
                "-Y",
                f"{Path(inputs.bwa_index.output.index.path) / index_fasta_name}",
                "input_reads_mate1.fastq.gz",
                "input_reads_mate2.fastq.gz",
            ]
            (Cmd["bwa-mem2"]["mem"][bwa_inputs] > aligned_sam)()
            self.progress(0.2)

        else:
            if inputs.aligned_reads.output.species != inputs.bwa_index.output.species:
                self.error(
                    "Species information for the input BAM file doesn't match the BWA-MEM2 index species information."
                )

            # Define output file names
            name = inputs.aligned_reads.entity_name
            if not name:
                bam_path = Path(inputs.aligned_reads.output.bam.path).name
                assert bam_path.endswith(".bam")
                name = bam_path[:-4]

            collate_inputs = [
                f"-@ {self.requirements.resources.cores}",
                "-O",
                inputs.aligned_reads.output.bam.path,
                "-",
            ]

            fastq_inputs = [
                f"-@ {self.requirements.resources.cores}",
                "-c 9",
                "-N",
                "-c singletons.fastq",
                "-",
            ]

            bwa_inputs = [
                "-K 100000000",
                "-v 3",
                f"-t {self.requirements.resources.cores}",
                "-p",
                "-Y",
                f"{Path(inputs.bwa_index.output.index.path) / index_fasta_name}",
                "-",
            ]

            (
                Cmd["samtools"]["collate"][collate_inputs]
                | Cmd["samtools"]["fastq"][fastq_inputs]
                | Cmd["bwa-mem2"]["mem"][bwa_inputs]
                > aligned_sam
            )()

        bam = f"{name}.bam"
        bam_stats = f"{name}_stats.txt"
        metrics_file = f"{name}_markduplicates_metrics.txt"

        # Convert aligned reads to BAM format
        # Samtools sort may require 4-5 GB RAM per thread, so the CPU
        # limit for this command is set to 4
        (
            Cmd["samtools"]["view"][
                "-1", "-@", min(4, self.requirements.resources.cores), aligned_sam
            ]
            > aligned_bam
        )()
        self.progress(0.25)

        # File cleanup
        Path(aligned_sam).unlink(missing_ok=True)

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
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["MarkDuplicates"][mark_duplicates_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("MarkDuplicates analysis failed.")

        self.progress(0.3)

        # File cleanup
        Path(aligned_bam).unlink(missing_ok=True)

        # Sort BAM file and fix NM and UQ tags
        sort_inputs = [
            "--INPUT",
            marked_dups,
            "--OUTPUT",
            sorted_temp,
            "--TMP_DIR",
            TMPDIR,
            "--SORT_ORDER",
            "coordinate",
            "--CREATE_INDEX",
            "false",
        ]

        return_code, _, _ = Cmd["gatk"]["SortSam"][sort_inputs] & TEE(retcode=None)
        if return_code:
            self.error("SortSam analysis failed.")

        self.progress(0.35)

        set_tag_inputs = [
            "--INPUT",
            sorted_temp,
            "--OUTPUT",
            sorted_bam,
            "--TMP_DIR",
            TMPDIR,
            "--CREATE_INDEX",
            "true",
            "--REFERENCE_SEQUENCE",
            inputs.ref_seq.output.fasta.path,
        ]

        return_code, _, _ = Cmd["gatk"]["SetNmMdAndUqTags"][set_tag_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("SetNmMdAndUqTags analysis failed.")

        self.progress(0.4)

        # File cleanup
        Path(marked_dups).unlink(missing_ok=True)
        Path(sorted_temp).unlink(missing_ok=True)

        # Set the read group information (required by BaseRecalibrator)
        rg_inputs = [
            "--INPUT",
            sorted_bam,
            "--VALIDATION_STRINGENCY",
            "STRICT",
            "--OUTPUT",
            sorted_rg,
            "--TMP_DIR",
            TMPDIR,
            "--RGLB",
            "WGS",
            "--RGPL",
            "ILLUMINA",
            "--RGPU",
            "X",
            "--RGSM",
            name,
            "--CREATE_INDEX",
            True,
        ]

        return_code, _, _ = Cmd["gatk"]["AddOrReplaceReadGroups"][rg_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("AddOrReplaceReadGroups tool failed.")

        self.progress(0.45)

        # File cleanup
        Path(sorted_bam).unlink(missing_ok=True)

        # Prepare files for scattering over chromosomes
        intervals_path = Path("intervals_folder")
        intervals_path.mkdir(exist_ok=True)

        if inputs.advanced_options.n_jobs:
            n_jobs = max(inputs.advanced_options.n_jobs, 1)
        else:
            n_jobs = max(self.requirements.resources.cores, 1)

        chromosome_sizes = "chromosome_sizes.bed"
        prepare_chromosome_sizes(
            fai_path=inputs.ref_seq.output.fai.path, bed_path=chromosome_sizes
        )

        split_intervals_inputs = [
            "-R",
            inputs.ref_seq.output.fasta.path,
            "-L",
            chromosome_sizes,
            "--scatter-count",
            n_jobs,
            "--subdivision-mode",
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
            "-O",
            str(intervals_path),
        ]
        return_code, _, _ = Cmd["gatk"]["SplitIntervals"][split_intervals_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("SplitIntervals tool failed.")

        # BaseRecalibrator
        recal_dir = Path("recal_tables")
        recal_dir.mkdir(exist_ok=True)

        intervals = [path for path in intervals_path.glob("*.interval_list")]
        return_codes = Parallel(n_jobs=n_jobs)(
            run_base_recalibration(
                intput_bam=sorted_rg,
                known_sites=inputs.known_sites,
                interval_path=interval_path,
                ref_seq_path=inputs.ref_seq.output.fasta.path,
                tmp=TMPDIR,
                parent_dir=recal_dir,
            )
            for interval_path in intervals
        )
        if any(return_codes):
            self.error("GATK BaseRecalibrator tool failed.")

        gather_bqsr_inputs = [
            "-O",
            recal_table,
            *prepare_scattered_inputs(results_dir=recal_dir, pattern="*.csv"),
        ]

        return_code, stdout, stderr = Cmd["gatk"]["GatherBQSRReports"][
            gather_bqsr_inputs
        ] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("GatherBQSRReports tool failed.")

        self.progress(0.6)

        # ApplyBQSR
        bqsr_dir = Path("bqsr_bams")
        bqsr_dir.mkdir(exist_ok=True)

        return_codes = Parallel(n_jobs=n_jobs)(
            run_apply_bqsr(
                intput_bam=sorted_rg,
                recal_table=recal_table,
                interval_path=interval_path,
                ref_seq_path=inputs.ref_seq.output.fasta.path,
                tmp=TMPDIR,
                parent_dir=bqsr_dir,
            )
            for interval_path in intervals
        )
        if any(return_codes):
            self.error("GATK ApplyBQSR tool failed.")

        gather_bam_inputs = [
            "-O",
            bam,
            *prepare_scattered_inputs(results_dir=bqsr_dir, pattern="*.bam"),
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, _, _ = Cmd["gatk"]["GatherBamFiles"][gather_bam_inputs] & TEE(
            retcode=None
        )
        if return_code:
            self.error("GatherBamFiles tool failed.")

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

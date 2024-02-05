"""Preprocess BAM file for RNA-seq variant calling."""

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


def prepare_read_groups(read_groups, error):
    """Prepare input read groups for GATK AddOrReplaceReadGroups."""
    input_groups = []
    present_tags = []
    for x in read_groups.split(";"):
        split_tag = x.split("=")
        input_groups.extend(split_tag)
        present_tags.append(split_tag[0])

    # Make sure all arguments to read_group are valid.
    all_tags = {
        "-LB",
        "-PL",
        "-PU",
        "-SM",
        "-CN",
        "-DS",
        "-DT",
        "-FO",
        "-ID",
        "-KS",
        "-PG",
        "-PI",
        "-PM",
        "-SO",
    }
    present_tag_set = set(present_tags)
    check_all_tags = present_tag_set.issubset(all_tags)
    if not check_all_tags:
        error(
            "One or more read_group argument(s) improperly formatted."
            f"{present_tag_set}"
        )

    # Check that there are no double entries of arguments to read_group.
    if len(present_tag_set) != len(present_tags):
        error("You have duplicate tags in read_group argument.")

    # Check that all mandatory arguments to read_group are present.
    mandatory_tags = {"-LB", "-PL", "-PU", "-SM"}
    check_tags = mandatory_tags.issubset(present_tag_set)
    if not check_tags:
        error(
            "Missing mandatory read_group argument(s) (-PL, -LB, -PU and -SM are mandatory)."
        )

    return input_groups


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
def run_split_ncigar_reads(
    input_bam, interval_path, ref_seq_path, tmp, parent_dir, memory
):
    """Run SplitNCigarReads on a specifed interval."""
    splitncigar_interval_bam = f"{parent_dir.name}/{interval_path.stem}.bam"

    apply_splitncigar_inputs = [
        "--java-options",
        f"-Xmx{memory}m",
        "-R",
        ref_seq_path,
        "-I",
        input_bam,
        "-O",
        splitncigar_interval_bam,
        "-L",
        interval_path,
        "--tmp-dir",
        tmp,
    ]
    return_code, stdout, stderr = Cmd["gatk"]["SplitNCigarReads"][
        apply_splitncigar_inputs
    ] & TEE(retcode=None)
    if return_code:
        print(f"Error in {interval_path.stem} interval.", stdout, stderr)
    return return_code


class RNASeqVC_Preprocess(Process):
    """Prepare BAM file from STAR aligner for HaplotypeCaller.

    This process includes steps MarkDuplicates, SplitNCigarReads,
    read-group assignment and base quality recalibration (BQSR).
    """

    slug = "rnaseq-vc-preprocess"
    name = "RNA-seq variant calling preprocess"
    category = "GATK"
    process_type = "data:alignment:bam:rnaseqvc"
    version = "1.3.0"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1"}
        },
        "resources": {
            "cores": 4,
            "memory": 65536,
            "storage": 600,
        },
    }
    entity = {"type": "sample"}

    data_name = "{{ bam|name|default('?') }}"

    class Input:
        """Input fields for RNASeqVC_Preprocess."""

        bam = DataField(
            data_type="alignment:bam",
            label="Alignment BAM file from STAR alignment",
        )
        ref_seq = DataField(
            data_type="seq:nucleotide",
            label="Reference sequence FASTA file",
        )
        known_sites = ListField(
            DataField(
                data_type="variants:vcf",
            ),
            label="List of known sites of variation",
            description="One or more databases of known polymorphic sites used to exclude regions around known "
            "polymorphisms from analysis.",
        )
        read_group = StringField(
            label="Replace read groups in BAM",
            description="Replace read groups in a BAM file. This argument enables the user to replace all read "
            "groups in the INPUT file with a single new read group and assign all reads to this read group in "
            "the OUTPUT BAM file. Addition or replacement is performed using GATK "
            "AddOrReplaceReadGroups tool. Input should take the form of -name=value delimited by a "
            '";", e.g. "-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1". See tool\'s '
            "documentation for more information on tag names. Note that PL, LB, PU and SM are require "
            "fields. See caveats of rewriting read groups in the documentation.",
            default="-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1",
        )

        class Advanced:
            """Advanced options."""

            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields for RNASeqVC_Preprocess."""

        bam = FileField(label="Preprocessed BAM file")
        bai = FileField(label="Index of BAM file")
        stats = FileField(label="Alignment statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")
        metrics_file = FileField(label="Metrics from MarkDuplicate process")

    def run(self, inputs, outputs):
        """Run analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        marked_dups = "marked_duplicates.bam"
        marked_dups_index = "marked_duplicates.bai"
        read_groups_file = "read_groups.bam"
        read_groups_index = "read_groups.bai"
        splitncigar = "splitNcigar.bam"
        recal_table = "recal_data.csv"

        file_name = Path(inputs.bam.output.bam.path).stem
        bam = f"{file_name}.bam"
        bai = f"{file_name}.bai"
        metrics_file = f"{file_name}_markduplicates_metrics.txt"

        read_groups = prepare_read_groups(
            read_groups=inputs.read_group, error=self.error
        )

        gc_threads = min(
            self.requirements.resources.cores, inputs.advanced.java_gc_threads
        )

        md_inputs = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--INPUT",
            inputs.bam.output.bam.path,
            "--VALIDATION_STRINGENCY",
            "STRICT",
            "--OUTPUT",
            marked_dups,
            "--METRICS_FILE",
            metrics_file,
            "--TMP_DIR",
            TMPDIR,
            "--CREATE_INDEX",
            "true",
        ]

        return_code, stdout, stderr = Cmd["gatk"]["MarkDuplicates"][md_inputs] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("MarkDuplicates analysis failed.")

        # Prepare files for scattering over chromosomes
        intervals_path = Path("intervals_folder")
        intervals_path.mkdir(exist_ok=True)

        n_jobs = max(inputs.advanced.java_gc_threads, self.requirements.resources.cores)

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
            intervals_path,
        ]
        return_code, stdout, stderr = Cmd["gatk"]["SplitIntervals"][
            split_intervals_inputs
        ] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("SplitIntervals tool failed.")

        # Prepare folder for SplitNCigarReads BAM files
        output_bams = Path("output_bams")
        output_bams.mkdir()

        memory = int(0.9 * (self.requirements.resources.memory / n_jobs))
        intervals = [path for path in intervals_path.glob("*.interval_list")]
        return_codes = Parallel(n_jobs=n_jobs)(
            run_split_ncigar_reads(
                input_bam=marked_dups,
                interval_path=interval_path,
                ref_seq_path=inputs.ref_seq.output.fasta.path,
                tmp=TMPDIR,
                parent_dir=output_bams,
                memory=memory,
            )
            for interval_path in intervals
        )
        if any(return_codes):
            self.error("GATK SplitNCigarReads tool failed.")

        Path(marked_dups).unlink(missing_ok=True)
        Path(marked_dups_index).unlink(missing_ok=True)

        input_lists = prepare_scattered_inputs(results_dir=output_bams, pattern="*.bam")
        gather_bam_inputs = [
            "-O",
            splitncigar,
            input_lists,
            "--TMP_DIR",
            TMPDIR,
        ]

        return_code, stdout, stderr = Cmd["gatk"]["GatherBamFiles"][
            gather_bam_inputs
        ] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("GatherBamFiles tool failed.")

        # Delete all files produces by SplitNCigarReads
        for bam_file in output_bams.glob("*"):
            Path(bam_file).unlink(missing_ok=True)
        Path(output_bams).rmdir()

        arg_rg = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--INPUT",
            splitncigar,
            "--VALIDATION_STRINGENCY",
            "STRICT",
            "--OUTPUT",
            read_groups_file,
            "--TMP_DIR",
            TMPDIR,
            "--CREATE_INDEX",
            "true",
        ]

        arg_rg.extend(read_groups)

        return_code, stdout, stderr = Cmd["gatk"]["AddOrReplaceReadGroups"][
            arg_rg
        ] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error("AddOrReplaceReadGroups failed.")

        # Delete merged BAM file produced by GatherBamFiles
        Path(splitncigar).unlink(missing_ok=True)

        recal_table = "recalibration.table"
        br_inputs = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--input",
            read_groups_file,
            "--output",
            recal_table,
            "--reference",
            inputs.ref_seq.output.fasta.path,
            "--read-validation-stringency",
            "STRICT",
            "--use-original-qualities",
            "--tmp-dir",
            TMPDIR,
        ]

        # Add known sites to the input parameters of BaseRecalibrator.
        for site in inputs.known_sites:
            br_inputs.extend(["--known-sites", f"{site.output.vcf.path}"])

        return_code, stdout, stderr = Cmd["gatk"]["BaseRecalibrator"][br_inputs] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("BaseRecalibrator failed.")

        # Apply base recalibration.
        ab_inputs = [
            "--java-options",
            f"-XX:ParallelGCThreads={gc_threads} -Xmx{inputs.advanced.max_heap_size}g",
            "--input",
            read_groups_file,
            "--output",
            bam,
            "--reference",
            inputs.ref_seq.output.fasta.path,
            "--bqsr-recal-file",
            recal_table,
            "--read-validation-stringency",
            "STRICT",
            "--use-original-qualities",
            "--tmp-dir",
            TMPDIR,
        ]

        return_code, stdout, stderr = Cmd["gatk"]["ApplyBQSR"][ab_inputs] & TEE(
            retcode=None
        )
        if return_code:
            print(stdout, stderr)
            self.error("ApplyBQSR failed.")

        # Delete BAM file produced by AddOrReplaceReadGroups
        Path(read_groups_file).unlink(missing_ok=True)
        Path(read_groups_index).unlink(missing_ok=True)

        stats = f"{bam}_stats.txt"
        (Cmd["samtools"]["flagstat"][bam] > stats)()

        outputs.bam = bam
        outputs.bai = bai
        outputs.stats = stats
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build
        outputs.metrics_file = metrics_file

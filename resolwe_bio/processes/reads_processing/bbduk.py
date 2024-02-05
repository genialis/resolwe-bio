"""Run BBDuk."""

import gzip
import shutil
from pathlib import Path

from joblib import Parallel, delayed, wrap_non_picklable_objects
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def prepare_inputs(mate1, mate2=None):
    """Prepare list with input, output and stats files."""
    input_reads = []

    if mate2:
        for fastq, fastq2 in zip(mate1, mate2):
            name, new_path, reads_basename = rename_input_files(fastq=fastq)
            name_mate2, new_path_mate2, reads_mate2_basename = rename_input_files(
                fastq=fastq2
            )

            lane_files = {
                "input_fastq": new_path,
                "input_fastq2": new_path_mate2,
                "output_fastq": f"{name}_preprocessed.fastq.gz",
                "output_fastq2": f"{name_mate2}_preprocessed.fastq.gz",
                "stats": f"{name}_statistics.txt",
                "original_name": f"{reads_basename}_preprocessed.fastq.gz",
                "original_name2": f"{reads_mate2_basename}_preprocessed.fastq.gz",
            }
            input_reads.append(lane_files)

    else:
        for fastq in mate1:
            name, new_path, reads_basename = rename_input_files(fastq=fastq)

            lane_files = {
                "input_fastq": new_path,
                "output_fastq": f"{name}_preprocessed.fastq.gz",
                "stats": f"{name}_statistics.txt",
                "original_name": f"{reads_basename}_preprocessed.fastq.gz",
            }
            input_reads.append(lane_files)

    return input_reads


def rename_input_files(fastq):
    """Rename input files."""
    reads_basename = Path(fastq.path)
    shutil.copy(reads_basename, Path.cwd())

    new_path = Path(reads_basename.name.replace(" ", "_"))
    Path(reads_basename.name).rename(new_path)
    assert new_path.name.endswith(".fastq.gz")
    name = new_path.name[:-9]

    return name, new_path, reads_basename.name[:-9]


def rename_preprocessed_files(input_files, paired_end=None):
    """Rename preprocessed files back to the original name."""
    for lane in input_files:
        Path(lane["output_fastq"]).rename(lane["original_name"])

        if paired_end:
            Path(lane["output_fastq2"]).rename(lane["original_name2"])


def prepare_fastqc(fastqgz, error):
    """Prepare FastQC data for output."""
    fastqc = []
    fastqc_url = []
    for fq in fastqgz:
        reads_name = Path(fq).name.replace(".fastq.gz", "")
        report_dir = Path("fastqc") / Path(f"{reads_name}_fastqc")

        fastqc_zip = Path(f"{reads_name}_fastqc.zip")
        if not fastqc_zip.is_file():
            error(f"FastQC failed to produce {fastqc_zip} file.")
        fastqc.append(str(fastqc_zip))

        fastqc_url.append(
            {
                "file": str(report_dir / "fastqc_report.html"),
                "refs": [str(report_dir)],
            }
        )
    return fastqc, fastqc_url


@delayed
@wrap_non_picklable_objects
def run_bbduk(input_reads, bbduk_inputs, paired_end=False):
    """Run BBDuk on seperate lanes."""

    if paired_end:
        input_file = [
            f"in='{input_reads['input_fastq']}'",
            f"in2='{input_reads['input_fastq2']}'",
            f"out='{input_reads['output_fastq']}'",
            f"out2='{input_reads['output_fastq2']}'",
            f"stats='{input_reads['stats']}'",
        ]
    else:
        input_file = [
            f"in='{input_reads['input_fastq']}'",
            f"out='{input_reads['output_fastq']}'",
            f"stats='{input_reads['stats']}'",
        ]

    bbduk_inputs = input_file + bbduk_inputs

    return_code, stdout, stderr = Cmd["bbduk.sh"][bbduk_inputs] & TEE(retcode=None)
    if return_code:
        print(stderr, stdout)
    return return_code, stderr


class BBDukSingle(Process):
    """Run BBDuk on single-end reads.

    BBDuk combines the most common data-quality-related trimming, filtering,
    and masking operations into a single high-performance tool. It is capable
    of quality-trimming and filtering, adapter-trimming, contaminant-filtering
    via kmer matching, sequence masking, GC-filtering, length filtering,
    entropy-filtering, format conversion, histogram generation, subsampling,
    quality-score recalibration, kmer cardinality estimation, and various
    other operations in a single pass. See
    [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
    for more information.
    """

    slug = "bbduk-single"
    name = "BBDuk (single-end)"
    process_type = "data:reads:fastq:single:bbduk"
    version = "3.1.2"
    category = "FASTQ processing"
    data_name = "{{ reads|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 8192,
        },
    }

    class Input:
        """Input fields to process BBDukSingle."""

        reads = DataField("reads:fastq:single", label="Reads")
        min_length = IntegerField(
            label="Minimum length",
            default=10,
            description="Reads shorter than the minimum length will be discarded after trimming.",
        )

        class Reference:
            """Reference."""

            sequences = ListField(
                DataField("seq:nucleotide"),
                label="Sequences",
                required=False,
                description="Reference sequences include adapters, contaminants, and "
                "degenerate sequences. They can be provided in a multi-sequence FASTA "
                "file or as a set of literal sequences below.",
            )
            literal_sequences = ListField(
                StringField(),
                label="Literal sequences",
                required=False,
                default=[],
                description="Literal sequences can be specified by inputting them one by "
                "one and pressing Enter after each sequence.",
            )

        class Processing:
            """Processing parameters."""

            kmer_length = IntegerField(
                label="Kmer length",
                default=27,
                description="Kmer length used for finding contaminants. "
                "Contaminants shorter than kmer length will not be found. "
                "Kmer length must be at least 1.",
            )
            check_reverse_complements = BooleanField(
                label="Check reverse complements",
                description="Look for reverse complements of kmers in addition to forward kmers",
                default=True,
            )
            mask_middle_base = BooleanField(
                label="Mask the middle base of a kmer",
                description="Treat the middle base of a kmer as a wildcard to increase "
                "sensitivity in the presence of errors.",
                default=True,
            )
            min_kmer_hits = IntegerField(
                label="Minimum number of kmer hits",
                default=1,
                description="Reads need at least this many matching kmers to be considered "
                "matching the reference.",
            )
            min_kmer_fraction = FloatField(
                label="Minimum kmer fraction",
                default=0.0,
                description="A read needs at least this fraction of its total kmers "
                "to hit a reference in order to be considered a match. If this and "
                "'Minimum number of kmer hits' are set, the greater is used.",
            )
            min_coverage_fraction = FloatField(
                label="Minimum coverage fraction",
                default=0.0,
                description="A read needs at least this fraction of its total bases to be "
                "covered by reference kmers to be considered a match. If specified, "
                "'Minimum coverage fraction' overrides 'Minimum number of kmer hits' and "
                "'Minimum kmer fraction'.",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for kmers (substitutions only)",
                description="Hamming distance i.e. the number of mismatches allowed in the kmer.",
                default=0,
            )
            query_hamming_distance = IntegerField(
                label="Hamming distance for query kmers",
                description="Set a hamming distance for query kmers instead of the read kmers. "
                "This makes the read processing much slower, but does not use additional memory.",
                default=0,
            )
            edit_distance = IntegerField(
                label="Maximum edit distance from reference kmers "
                "(substitutions and indels)",
                default=0,
            )
            hamming_distance2 = IntegerField(
                label="Hamming distance for short kmers when looking for shorter kmers",
                default=0,
            )
            query_hamming_distance2 = IntegerField(
                label="Hamming distance for short query kmers when looking for shorter kmers",
                default=0,
            )
            edit_distance2 = IntegerField(
                label="Maximum edit distance from short reference kmers "
                "(substitutions and indels) when looking for shorter kmers",
                default=0,
            )
            forbid_N = BooleanField(
                label="Forbid matching of read kmers containing N",
                default=False,
                description="By default, these will match a reference 'A' if"
                "'Maximum Hamming distance for kmers' > 0 or "
                "'Maximum edit distance from reference kmers' > 0, to increase sensitivity.",
            )
            find_best_match = BooleanField(
                label="Find best match",
                description="If multiple matches, associate read with sequence sharing most kmers.",
                default=True,
            )

        class Operations:
            """Trimming, filtering and masking parameters."""

            k_trim = StringField(
                label="Trimming protocol to remove bases matching reference kmers from reads",
                choices=[
                    ("f", "Don't trim"),
                    ("r", "Trim to the right"),
                    ("l", "Trim to the left"),
                ],
                default="f",
            )
            k_mask = StringField(
                label="Symbol to replace bases matching reference kmers",
                default="f",
                description="Allows any non-whitespace character other than t or f. "
                "Processes short kmers on both ends.",
            )
            mask_fully_covered = BooleanField(
                label="Only mask bases that are fully covered by kmers",
                default=False,
            )
            min_k = IntegerField(
                label="Look for shorter kmers at read tips down to this length when k-trimming "
                "or masking",
                default=-1,
                description="-1 means disabled. Enabling this will disable treating the middle "
                "base of a kmer as a wildcard to increase sensitivity in the presence of errors.",
            )
            quality_trim = StringField(
                label="Trimming protocol to remove bases with quality below "
                "the minimum average region quality from read ends",
                choices=[
                    ("f", "Trim neither end"),
                    ("rl", "Trim both ends"),
                    ("r", "Trim only right end"),
                    ("l", "Trim only left end"),
                    ("w", "Use sliding window"),
                ],
                default="f",
                description="Performed after looking for kmers. If enabled, set also "
                "'Average quality below which to trim region'.",
            )
            trim_quality = IntegerField(
                label="Average quality below which to trim region",
                default=6,
                disabled="operations.quality_trim === 'f'",
                description="Set trimming protocol to enable this parameter.",
            )
            quality_encoding_offset = StringField(
                label="Quality encoding offset",
                choices=[
                    ("33", "Sanger / Illumina 1.8+ (33)"),
                    ("64", "Illumina up to 1.3+, 1.5+ (64)"),
                    ("auto", "Auto"),
                ],
                default="auto",
                description="Quality encoding offset for input FASTQ files.",
            )
            ignore_bad_quality = BooleanField(
                label="Don't crash if quality values appear to be incorrect",
                default=False,
            )
            trim_poly_A = IntegerField(
                label="Minimum length of poly-A or poly-T tails to trim on either end of reads",
                default=0,
            )
            min_length_fraction = FloatField(
                label="Minimum length fraction",
                default=0.0,
                description="Reads shorter than this fraction of original length after "
                "trimming will be discarded.",
            )
            max_length = IntegerField(
                label="Maximum length",
                required=False,
                description="Reads longer than this after trimming will be discarded.",
            )
            min_average_quality = IntegerField(
                label="Minimum average quality",
                default=0,
                description="Reads with average quality (after trimming) below this will be discarded.",
            )
            min_average_quality_bases = IntegerField(
                label="Number of initial bases to calculate minimum average quality from",
                default=0,
                description="If positive, calculate minimum average quality "
                "from this many initial bases",
            )
            min_base_quality = IntegerField(
                label="Minimum base quality below which reads are discarded after trimming",
                default=0,
            )
            min_consecutive_bases = IntegerField(
                label="Minimum number of consecutive called bases",
                default=0,
            )
            trim_pad = IntegerField(
                label="Number of bases to trim around matching kmers",
                default=0,
            )
            min_overlap = IntegerField(
                label="Minum number of overlapping bases",
                default=14,
                description="Require this many bases of overlap for detection.",
            )
            min_insert = IntegerField(
                label="Minimum insert size",
                default=40,
                description="Require insert size of at least this for overlap. "
                "Should be reduced to 16 for small RNA sequencing.",
            )
            force_trim_left = IntegerField(
                label="Position from which to trim bases to the left",
                default=0,
            )
            force_trim_right = IntegerField(
                label="Position from which to trim bases to the right",
                default=0,
            )
            force_trim_right2 = IntegerField(
                label="Number of bases to trim from the right end",
                default=0,
            )
            force_trim_mod = IntegerField(
                label="Modulo to right-trim reads",
                default=0,
                description="Trim reads to the largest multiple of modulo.",
            )
            restrict_left = IntegerField(
                label="Number of leftmost bases to look in for kmer matches",
                default=0,
            )
            restrict_right = IntegerField(
                label="Number of rightmost bases to look in for kmer matches",
                default=0,
            )
            min_GC = FloatField(
                label="Minimum GC content",
                default=0.0,
                description="Discard reads with lower GC content.",
            )
            max_GC = FloatField(
                label="Maximum GC content",
                default=1.0,
                description="Discard reads with higher GC content.",
            )
            maxns = IntegerField(
                label="Max Ns after trimming",
                default=-1,
                description="If non-negative, reads with more Ns than this "
                "(after trimming) will be discarded.",
            )
            toss_junk = BooleanField(
                label="Discard reads with invalid characters as bases",
                default=False,
            )

        class HeaderParsing:
            """Header-parsing parameters."""

            chastity_filter = BooleanField(
                label="Discard reads that fail Illumina chastity filtering",
                default=False,
                description="Discard reads with id containing ' 1:Y:' or ' 2:Y:'.",
            )
            barcode_filter = BooleanField(
                label="Remove reads with unexpected barcodes",
                default=False,
                description="Remove reads with unexpected barcodes if barcodes are set, "
                "or barcodes containing 'N' otherwise. A barcode must be the last part "
                "of the read header.",
            )
            barcode_files = ListField(
                DataField("seq:nucleotide"),
                label="Barcode sequences",
                description="FASTA file(s) with barcode sequences.",
                required=False,
            )
            barcode_sequences = ListField(
                StringField(),
                label="Literal barcode sequences",
                required=False,
                default=[],
                description="Literal barcode sequences can be specified by inputting "
                "them one by one and pressing Enter after each sequence.",
            )
            x_min = IntegerField(
                label="Minimum X coordinate",
                default=-1,
                description="If positive, discard reads with a smaller X coordinate.",
            )
            y_min = IntegerField(
                label="Minimum Y coordinate",
                default=-1,
                description="If positive, discard reads with a smaller Y coordinate.",
            )
            x_max = IntegerField(
                label="Maximum X coordinate",
                default=-1,
                description="If positive, discard reads with a larger X coordinate.",
            )
            y_max = IntegerField(
                label="Maximum Y coordinate",
                default=-1,
                description="If positive, discard reads with a larger Y coordinate.",
            )

        class Complexity:
            """Complexity parameters."""

            entropy = FloatField(
                label="Minimum entropy",
                default=-1.0,
                description="Set between 0 and 1 to filter reads with entropy below that value. "
                "Higher is more stringent.",
            )
            entropy_window = IntegerField(
                label="Length of sliding window used to calculate entropy",
                default=50,
                description="To use the sliding window set minimum entropy in range between 0.0 and 1.0.",
            )
            entropy_k = IntegerField(
                label="Length of kmers used to calcuate entropy",
                default=5,
            )
            entropy_mask = BooleanField(
                label="Mask low-entropy parts of sequences with N instead of discarding",
                default=False,
            )
            min_base_frequency = IntegerField(
                label="Minimum base frequency",
                default=0,
            )

        class Fastqc:
            """FastQC parameters."""

            nogroup = BooleanField(
                label="Disable grouping of bases for reads >50bp",
                default=False,
                description="All reports will show data for every base in the read. Using this option "
                "will cause fastqc to crash and burn if you use it on really long reads.",
            )

        reference = GroupField(Reference, label="Reference")

        processing = GroupField(Processing, label="Processing parameters")

        operations = GroupField(
            Operations, label="Trimming, filtering and masking parameters."
        )

        header_parsing = GroupField(HeaderParsing, label="Header-parsing parameters")

        complexity = GroupField(Complexity, label="Complexity parameters")

        fastqc = GroupField(Fastqc, label="FastQC parameters")

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Remaining reads")
        statistics = ListField(FileField(), label="Statistics")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC")
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run analysis."""

        input_references = "input_references.fasta.gz"
        input_barcodes = "input_barcodes.fasta.gz"

        num_of_lanes = len(inputs.reads.output.fastq)
        input_reads = prepare_inputs(mate1=inputs.reads.output.fastq)

        if inputs.reference.sequences:
            with gzip.open(input_references, "wb") as outfile:
                for sequences in inputs.reference.sequences:
                    with gzip.open(sequences.output.fastagz.path, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)

        if inputs.header_parsing.barcode_files:
            with gzip.open(input_barcodes, "wb") as outfile:
                for barcode_file in inputs.header_parsing.barcode_files:
                    with gzip.open(barcode_file.output.fastagz.path, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            barcodes = [input_barcodes] + inputs.header_parsing.barcode_sequences
        else:
            barcodes = inputs.header_parsing.barcode_sequences

        self.progress(0.1)

        args = [
            "statscolumns=5",
            f"k={inputs.processing.kmer_length}",
            f"rcomp={inputs.processing.check_reverse_complements}",
            f"maskmiddle={inputs.processing.mask_middle_base}",
            f"minkmerhits={inputs.processing.min_kmer_hits}",
            f"minkmerfraction={inputs.processing.min_kmer_fraction}",
            f"mincovfraction={inputs.processing.min_coverage_fraction}",
            f"hammingdistance={inputs.processing.hamming_distance}",
            f"qhdist={inputs.processing.query_hamming_distance}",
            f"editdistance={inputs.processing.edit_distance}",
            f"hammingdistance2={inputs.processing.hamming_distance2}",
            f"qhdist2={inputs.processing.query_hamming_distance2}",
            f"editdistance2={inputs.processing.edit_distance2}",
            f"forbidn={inputs.processing.forbid_N}",
            f"findbestmatch={inputs.processing.find_best_match}",
            f"maskfullycovered={inputs.operations.mask_fully_covered}",
            f"mink={inputs.operations.min_k}",
            f"trimq={inputs.operations.trim_quality}",
            f"qin={inputs.operations.quality_encoding_offset}",
            f"trimpolya={inputs.operations.trim_poly_A}",
            f"minlength={inputs.min_length}",
            f"minlengthfraction={inputs.operations.min_length_fraction}",
            f"minavgquality={inputs.operations.min_average_quality}",
            f"maqb={inputs.operations.min_average_quality_bases}",
            f"minbasequality={inputs.operations.min_base_quality}",
            f"minconsecutivebases={inputs.operations.min_consecutive_bases}",
            f"trimpad={inputs.operations.trim_pad}",
            f"minoverlap={inputs.operations.min_overlap}",
            f"mininsert={inputs.operations.min_insert}",
            f"forcetrimleft={inputs.operations.force_trim_left}",
            f"forcetrimright={inputs.operations.force_trim_right}",
            f"forcetrimright2={inputs.operations.force_trim_right2}",
            f"forcetrimmod={inputs.operations.force_trim_mod}",
            f"restrictleft={inputs.operations.restrict_left}",
            f"restrictright={inputs.operations.restrict_right}",
            f"mingc={inputs.operations.min_GC}",
            f"maxgc={inputs.operations.max_GC}",
            f"maxns={inputs.operations.maxns}",
            f"tossjunk={inputs.operations.toss_junk}",
            f"chastityfilter={inputs.header_parsing.chastity_filter}",
            f"barcodefilter={inputs.header_parsing.barcode_filter}",
            f"xmin={inputs.header_parsing.x_min}",
            f"ymin={inputs.header_parsing.y_min}",
            f"xmax={inputs.header_parsing.x_max}",
            f"ymax={inputs.header_parsing.y_max}",
            f"entropy={inputs.complexity.entropy}",
            f"entropywindow={inputs.complexity.entropy_window}",
            f"entropyk={inputs.complexity.entropy_k}",
            f"minbasefrequency={inputs.complexity.min_base_frequency}",
            f"entropymask={inputs.complexity.entropy_mask}",
            f"-Xmx{int(0.85*(self.requirements.resources.memory/num_of_lanes))}m",
        ]

        if self.requirements.resources.cores >= num_of_lanes:
            args.append(
                f"threads={int(self.requirements.resources.cores//num_of_lanes)}"
            )
            n_jobs = num_of_lanes
        else:
            self.warning(
                f"There are more sequencing lanes ({num_of_lanes}) than there are "
                f"available cores ({self.requirements.resources.cores}). For the "
                "most optimal performance, use at least the same number of lanes "
                "and cores."
            )
            args.append("threads=1")
            n_jobs = self.requirements.resources.cores

        if inputs.reference.sequences:
            args.append(f"ref={input_references}")

        if inputs.reference.literal_sequences:
            literal_sequences_joined = ",".join(inputs.reference.literal_sequences)
            args.append(f"literal={literal_sequences_joined}")

        if inputs.operations.k_trim != "f":
            args.append(f"ktrim={inputs.operations.k_trim}")

        if inputs.operations.k_mask != "f" and inputs.operations.k_mask != "t":
            args.append(f"kmask={inputs.operations.k_mask}")

        if inputs.operations.quality_trim != "f":
            args.append(f"qtrim={inputs.operations.quality_trim}")

        if inputs.operations.ignore_bad_quality:
            args.append("ignorebadquality")

        if inputs.operations.max_length:
            args.append(f"maxlength={inputs.operations.max_length}")

        if barcodes:
            barcodes = ",".join(barcodes)
            args.append(f"barcodes={barcodes}")

        process_outputs = Parallel(n_jobs=n_jobs)(
            run_bbduk(input_reads=input_set, bbduk_inputs=args)
            for input_set in input_reads
        )

        for output in process_outputs:
            if output[0]:
                self.error("BBDuk failed.", output[1])

        self.progress(0.7)

        statistics = []
        for stats in input_reads:
            with open(stats["stats"], "rb") as orig_file:
                with gzip.open(f"{stats['stats']}.gz", "wb") as zipped_file:
                    zipped_file.writelines(orig_file)
            statistics.append(f"{stats['stats']}.gz")

        output_path = Path("./fastqc")
        output_path.mkdir(exist_ok=True)

        rename_preprocessed_files(input_files=input_reads)

        fastqgz = [fastq["original_name"] for fastq in input_reads]
        fastqc_inputs = fastqgz + ["--extract", f"--outdir={str(output_path)}"]

        if inputs.fastqc.nogroup:
            fastqc_inputs.append("--no-group")

        return_code, _, stderr = Cmd["fastqc"][fastqc_inputs] & TEE(retcode=None)
        if return_code:
            self.error("FastQC failed. ", stderr)

        for fastqc_zip in output_path.glob("*_fastqc.zip"):
            shutil.move(str(fastqc_zip), ".")

        fastqc, fastqc_url = prepare_fastqc(fastqgz=fastqgz, error=self.error)

        if (
            inputs.operations.k_trim == "f"
            and inputs.operations.quality_trim == "f"
            and not inputs.reference.sequences
            and not inputs.reference.literal_sequences
            and inputs.operations.min_average_quality <= 0
        ):
            self.warning(
                "Reference sequences, trimming mode, and minimum average quality are "
                "unspecified. Only filtering of reads by length is applied."
            )

        outputs.fastq = fastqgz
        outputs.statistics = statistics
        outputs.fastqc_url = fastqc_url
        outputs.fastqc_archive = fastqc


class BBDukPaired(Process):
    """Run BBDuk on paired-end reads.

    BBDuk combines the most common data-quality-related trimming, filtering,
    and masking operations into a single high-performance tool. It is capable
    of quality-trimming and filtering, adapter-trimming, contaminant-filtering
    via kmer matching, sequence masking, GC-filtering, length filtering,
    entropy-filtering, format conversion, histogram generation, subsampling,
    quality-score recalibration, kmer cardinality estimation, and various
    other operations in a single pass. See
    [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
    for more information.
    """

    slug = "bbduk-paired"
    name = "BBDuk (paired-end)"
    process_type = "data:reads:fastq:paired:bbduk"
    version = "3.1.2"
    category = "FASTQ processing"
    data_name = "{{ reads|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 8192,
        },
    }

    class Input:
        """Input fields to process BBDukPaired."""

        reads = DataField("reads:fastq:paired", label="Reads")
        min_length = IntegerField(
            label="Minimum length",
            default=10,
            description="Reads shorter than the minimum length will be discarded after trimming.",
        )

        class Reference:
            """Reference."""

            sequences = ListField(
                DataField("seq:nucleotide"),
                label="Sequences",
                required=False,
                description="Reference sequences include adapters, contaminants, and "
                "degenerate sequences. They can be provided in a multi-sequence FASTA "
                "file or as a set of literal sequences below.",
            )
            literal_sequences = ListField(
                StringField(),
                label="Literal sequences",
                required=False,
                default=[],
                description="Literal sequences can be specified by inputting them one by "
                "one and pressing Enter after each sequence.",
            )

        class Processing:
            """Processing parameters."""

            kmer_length = IntegerField(
                label="Kmer length",
                default=27,
                description="Kmer length used for finding contaminants. "
                "Contaminants shorter than kmer length will not be found. "
                "Kmer length must be at least 1.",
            )
            check_reverse_complements = BooleanField(
                label="Check reverse complements",
                description="Look for reverse complements of kmers in addition to forward kmers.",
                default=True,
            )
            mask_middle_base = BooleanField(
                label="Mask the middle base of a kmer",
                description="Treat the middle base of a kmer as a wildcard to increase sensitivity "
                "in the presence of errors.",
                default=True,
            )
            min_kmer_hits = IntegerField(
                label="Minimum number of kmer hits",
                default=1,
                description="Reads need at least this many matching kmers to be considered "
                "matching the reference.",
            )
            min_kmer_fraction = FloatField(
                label="Minimum kmer fraction",
                default=0.0,
                description="A read needs at least this fraction of its total kmers "
                "to hit a reference in order to be considered a match. If this and "
                "'Minimum number of kmer hits' are set, the greater is used.",
            )
            min_coverage_fraction = FloatField(
                label="Minimum kmer fraction",
                default=0.0,
                description="A read needs at least this fraction of its total bases to be "
                "covered by reference kmers to be considered a match. If specified, "
                "'Minimum coverage fraction' overrides 'Minimum number of kmer hits' and "
                "'Minimum kmer fraction'.",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for kmers (substitutions only)",
                description="Hamming distance i.e. the number of mismatches allowed in the kmer.",
                default=0,
            )
            query_hamming_distance = IntegerField(
                label="Hamming distance for query kmers",
                description="Set a hamming distance for query kmers instead of the read kmers. "
                "This makes the read processing much slower, but does not use additional memory.",
                default=0,
            )
            edit_distance = IntegerField(
                label="Maximum edit distance from reference kmers "
                "(substitutions and indels)",
                default=0,
            )
            hamming_distance2 = IntegerField(
                label="Hamming distance for short kmers when looking for shorter kmers",
                default=0,
            )
            query_hamming_distance2 = IntegerField(
                label="Hamming distance for short query kmers when looking for shorter kmers",
                default=0,
            )
            edit_distance2 = IntegerField(
                label="Maximum edit distance from short reference kmers "
                "(substitutions and indels) when looking for shorter kmers",
                default=0,
            )
            forbid_N = BooleanField(
                label="Forbid matching of read kmers containing N",
                default=False,
                description="By default, these will match a reference 'A' if"
                "'Maximum Hamming distance for kmers' > 0 or "
                "'Maximum edit distance from reference kmers' > 0, to increase sensitivity.",
            )
            find_best_match = BooleanField(
                label="Find best match",
                description="If multiple matches, associate read with sequence sharing most kmers.",
                default=True,
            )
            remove_if_either_bad = BooleanField(
                label="Remove both sequences of a paired-end read, if either of them is to be removed",
                default=True,
            )
            perform_error_correction = BooleanField(
                label="Perform error correction with BBMerge prior to kmer operations",
                default=False,
            )

        class Operations:
            """Trimming, filtering and masking parameters."""

            k_trim = StringField(
                label="Trimming protocol to remove bases matching reference kmers from reads",
                choices=[
                    ("f", "Don't trim"),
                    ("r", "Trim to the right"),
                    ("l", "Trim to the left"),
                ],
                default="f",
            )
            k_mask = StringField(
                label="Symbol to replace bases matching reference kmers",
                default="f",
                description="Allows any non-whitespace character other than t or f. "
                "Processes short kmers on both ends.",
            )
            mask_fully_covered = BooleanField(
                label="Only mask bases that are fully covered by kmers",
                default=False,
            )
            min_k = IntegerField(
                label="Look for shorter kmers at read tips down to this length when k-trimming "
                "or masking",
                default=-1,
                description="-1 means disabled. Enabling this will disable treating the middle "
                "base of a kmer as a wildcard to increase sensitivity in the presence of errors.",
            )
            quality_trim = StringField(
                label="Trimming protocol to remove bases with quality below "
                "the minimum average region quality from read ends",
                choices=[
                    ("f", "Trim neither end"),
                    ("rl", "Trim both ends"),
                    ("r", "Trim only right end"),
                    ("l", "Trim only left end"),
                    ("w", "Use sliding window"),
                ],
                default="f",
                description="Performed after looking for kmers. If enabled, set also "
                "'Average quality below which to trim region'.",
            )
            trim_quality = IntegerField(
                label="Average quality below which to trim region",
                default=6,
                disabled="operations.quality_trim === 'f'",
                description="Set trimming protocol to enable this parameter.",
            )
            quality_encoding_offset = StringField(
                label="Quality encoding offset",
                choices=[
                    ("33", "Sanger / Illumina 1.8+ (33)"),
                    ("64", "Illumina up to 1.3+, 1.5+ (64)"),
                    ("auto", "Auto"),
                ],
                default="auto",
                description="Quality encoding offset for input FASTQ files.",
            )
            ignore_bad_quality = BooleanField(
                label="Don't crash if quality values appear to be incorrect",
                default=False,
            )
            trim_poly_A = IntegerField(
                label="Minimum length of poly-A or poly-T tails to trim on either end of reads",
                default=0,
            )
            min_length_fraction = FloatField(
                label="Minimum length fraction",
                default=0.0,
                description="Reads shorter than this fraction of original length after "
                "trimming will be discarded.",
            )
            max_length = IntegerField(
                label="Maximum length",
                required=False,
                description="Reads longer than this after trimming will be discarded.",
            )
            min_average_quality = IntegerField(
                label="Minimum average quality",
                default=0,
                description="Reads with average quality (after trimming) below this will be discarded.",
            )
            min_average_quality_bases = IntegerField(
                label="Number of initial bases to calculate minimum average quality from",
                default=0,
                description="If positive, calculate minimum average quality "
                "from this many initial bases",
            )
            min_base_quality = IntegerField(
                label="Minimum base quality below which reads are discarded after trimming",
                default=0,
            )
            min_consecutive_bases = IntegerField(
                label="Minimum number of consecutive called bases",
                default=0,
            )
            trim_pad = IntegerField(
                label="Number of bases to trim around matching kmers",
                default=0,
            )
            trim_by_overlap = BooleanField(
                label="Trim adapters based on where paired-end reads overlap",
                default=False,
            )
            strict_overlap = BooleanField(
                label="Adjust sensitivity in "
                "'Trim adapters based on where paired-end reads overlap' mode",
                default=True,
            )
            min_overlap = IntegerField(
                label="Minum number of overlapping bases",
                default=14,
                description="Require this many bases of overlap for detection.",
            )
            min_insert = IntegerField(
                label="Minimum insert size",
                default=40,
                description="Require insert size of at least this for overlap. "
                "Should be reduced to 16 for small RNA sequencing.",
            )
            trim_pairs_evenly = BooleanField(
                label="Trim both sequences of paired-end reads to the minimum "
                "length of either sequence",
                default=False,
            )
            force_trim_left = IntegerField(
                label="Position from which to trim bases to the left",
                default=0,
            )
            force_trim_right = IntegerField(
                label="Position from which to trim bases to the right",
                default=0,
            )
            force_trim_right2 = IntegerField(
                label="Number of bases to trim from the right end",
                default=0,
            )
            force_trim_mod = IntegerField(
                label="Modulo to right-trim reads",
                default=0,
                description="Trim reads to the largest multiple of modulo.",
            )
            restrict_left = IntegerField(
                label="Number of leftmost bases to look in for kmer matches",
                default=0,
            )
            restrict_right = IntegerField(
                label="Number of rightmost bases to look in for kmer matches",
                default=0,
            )
            min_GC = FloatField(
                label="Minimum GC content",
                default=0.0,
                description="Discard reads with lower GC content.",
            )
            max_GC = FloatField(
                label="Maximum GC content",
                default=1.0,
                description="Discard reads with higher GC content.",
            )
            maxns = IntegerField(
                label="Max Ns after trimming",
                default=-1,
                description="If non-negative, reads with more Ns than this "
                "(after trimming) will be discarded.",
            )
            toss_junk = BooleanField(
                label="Discard reads with invalid characters as bases",
                default=False,
            )

        class HeaderParsing:
            """Header-parsing parameters."""

            chastity_filter = BooleanField(
                label="Discard reads that fail Illumina chastity filtering",
                default=False,
                description="Discard reads with id containing ' 1:Y:' or ' 2:Y:'.",
            )
            barcode_filter = BooleanField(
                label="Remove reads with unexpected barcodes",
                default=False,
                description="Remove reads with unexpected barcodes if barcodes are set, "
                "or barcodes containing 'N' otherwise. A barcode must be the last part "
                "of the read header.",
            )
            barcode_files = ListField(
                DataField("seq:nucleotide"),
                label="Barcode sequences",
                description="FASTA file(s) with barcode sequences.",
                required=False,
            )
            barcode_sequences = ListField(
                StringField(),
                label="Literal barcode sequences",
                required=False,
                default=[],
                description="Literal barcode sequences can be specified by inputting "
                "them one by one and pressing Enter after each sequence.",
            )
            x_min = IntegerField(
                label="Minimum X coordinate",
                default=-1,
                description="If positive, discard reads with a smaller X coordinate.",
            )
            y_min = IntegerField(
                label="Minimum Y coordinate",
                default=-1,
                description="If positive, discard reads with a smaller Y coordinate.",
            )
            x_max = IntegerField(
                label="Maximum X coordinate",
                default=-1,
                description="If positive, discard reads with a larger X coordinate.",
            )
            y_max = IntegerField(
                label="Maximum Y coordinate",
                default=-1,
                description="If positive, discard reads with a larger Y coordinate.",
            )

        class Complexity:
            """Complexity parameters."""

            entropy = FloatField(
                label="Minimum entropy",
                default=-1.0,
                description="Set between 0 and 1 to filter reads with entropy below that value. "
                "Higher is more stringent.",
            )
            entropy_window = IntegerField(
                label="Length of sliding window used to calculate entropy",
                default=50,
                description="To use the sliding window set minimum entropy in range between 0.0 and 1.0.",
            )
            entropy_k = IntegerField(
                label="Length of kmers used to calcuate entropy",
                default=5,
            )
            entropy_mask = BooleanField(
                label="Mask low-entropy parts of sequences with N instead of discarding",
                default=False,
            )
            min_base_frequency = IntegerField(
                label="Minimum base frequency",
                default=0,
            )

        class Fastqc:
            """FastQC parameters."""

            nogroup = BooleanField(
                label="Disable grouping of bases for reads >50bp",
                default=False,
                description="All reports will show data for every base in the read. Using this option "
                "will cause fastqc to crash and burn if you use it on really long reads.",
            )

        reference = GroupField(Reference, label="Reference")

        processing = GroupField(Processing, label="Processing parameters")

        operations = GroupField(
            Operations, label="Trimming, filtering and masking parameters."
        )

        header_parsing = GroupField(HeaderParsing, label="Header-parsing parameters")

        complexity = GroupField(Complexity, label="Complexity parameters")

        fastqc = GroupField(Fastqc, label="FastQC parameters")

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Remaining upstream reads")
        fastq2 = ListField(FileField(), label="Remaining downstream reads")
        statistics = ListField(FileField(), label="Statistics")
        fastqc_url = ListField(
            FileHtmlField(), label="Upstream quality control with FastQC"
        )
        fastqc_url2 = ListField(
            FileHtmlField(), label="Downstream quality control with FastQC"
        )
        fastqc_archive = ListField(
            FileField(), label="Download upstream FastQC archive"
        )
        fastqc_archive2 = ListField(
            FileField(), label="Download downstream FastQC archive"
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        input_references = "input_references.fasta.gz"
        input_barcodes = "input_barcodes.fasta.gz"

        num_of_lanes = len(inputs.reads.output.fastq)
        input_reads = prepare_inputs(
            mate1=inputs.reads.output.fastq, mate2=inputs.reads.output.fastq2
        )

        if inputs.reference.sequences:
            with gzip.open(input_references, "wb") as outfile:
                for sequences in inputs.reference.sequences:
                    with gzip.open(sequences.output.fastagz.path, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)

        if inputs.header_parsing.barcode_files:
            with gzip.open(input_barcodes, "wb") as outfile:
                for barcode_file in inputs.header_parsing.barcode_files:
                    with gzip.open(barcode_file.output.fastagz.path, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
            barcodes = [input_barcodes] + inputs.header_parsing.barcode_sequences
        else:
            barcodes = inputs.header_parsing.barcode_sequences

        self.progress(0.1)

        args = [
            "statscolumns=5",
            f"k={inputs.processing.kmer_length}",
            f"rcomp={inputs.processing.check_reverse_complements}",
            f"maskmiddle={inputs.processing.mask_middle_base}",
            f"minkmerhits={inputs.processing.min_kmer_hits}",
            f"minkmerfraction={inputs.processing.min_kmer_fraction}",
            f"mincovfraction={inputs.processing.min_coverage_fraction}",
            f"hammingdistance={inputs.processing.hamming_distance}",
            f"qhdist={inputs.processing.query_hamming_distance}",
            f"editdistance={inputs.processing.edit_distance}",
            f"hammingdistance2={inputs.processing.hamming_distance2}",
            f"qhdist2={inputs.processing.query_hamming_distance2}",
            f"editdistance2={inputs.processing.edit_distance2}",
            f"forbidn={inputs.processing.forbid_N}",
            f"removeifeitherbad={inputs.processing.remove_if_either_bad}",
            f"findbestmatch={inputs.processing.find_best_match}",
            f"ecco={inputs.processing.perform_error_correction}",
            f"maskfullycovered={inputs.operations.mask_fully_covered}",
            f"mink={inputs.operations.min_k}",
            f"trimq={inputs.operations.trim_quality}",
            f"qin={inputs.operations.quality_encoding_offset}",
            f"trimpolya={inputs.operations.trim_poly_A}",
            f"minlength={inputs.min_length}",
            f"minlengthfraction={inputs.operations.min_length_fraction}",
            f"minavgquality={inputs.operations.min_average_quality}",
            f"maqb={inputs.operations.min_average_quality_bases}",
            f"minbasequality={inputs.operations.min_base_quality}",
            f"minconsecutivebases={inputs.operations.min_consecutive_bases}",
            f"trimpad={inputs.operations.trim_pad}",
            f"trimbyoverlap={inputs.operations.trim_by_overlap}",
            f"trimpairsevenly={inputs.operations.trim_pairs_evenly}",
            f"minoverlap={inputs.operations.min_overlap}",
            f"mininsert={inputs.operations.min_insert}",
            f"forcetrimleft={inputs.operations.force_trim_left}",
            f"forcetrimright={inputs.operations.force_trim_right}",
            f"forcetrimright2={inputs.operations.force_trim_right2}",
            f"forcetrimmod={inputs.operations.force_trim_mod}",
            f"restrictleft={inputs.operations.restrict_left}",
            f"restrictright={inputs.operations.restrict_right}",
            f"mingc={inputs.operations.min_GC}",
            f"maxgc={inputs.operations.max_GC}",
            f"maxns={inputs.operations.maxns}",
            f"tossjunk={inputs.operations.toss_junk}",
            f"chastityfilter={inputs.header_parsing.chastity_filter}",
            f"barcodefilter={inputs.header_parsing.barcode_filter}",
            f"xmin={inputs.header_parsing.x_min}",
            f"ymin={inputs.header_parsing.y_min}",
            f"xmax={inputs.header_parsing.x_max}",
            f"ymax={inputs.header_parsing.y_max}",
            f"entropy={inputs.complexity.entropy}",
            f"entropywindow={inputs.complexity.entropy_window}",
            f"entropyk={inputs.complexity.entropy_k}",
            f"minbasefrequency={inputs.complexity.min_base_frequency}",
            f"entropymask={inputs.complexity.entropy_mask}",
            f"-Xmx{int(0.85*(self.requirements.resources.memory/num_of_lanes))}m",
        ]

        if self.requirements.resources.cores >= num_of_lanes:
            args.append(
                f"threads={int(self.requirements.resources.cores//num_of_lanes)}"
            )
            n_jobs = num_of_lanes
        else:
            self.warning(
                f"There are more sequencing lanes ({num_of_lanes}) than there are "
                f"available cores ({self.requirements.resources.cores}). For the "
                "most optimal performance, use at least the same number of lanes "
                "and cores."
            )
            args.append("threads=1")
            n_jobs = self.requirements.resources.cores

        if inputs.reference.sequences:
            args.append(f"ref={input_references}")

        if inputs.reference.literal_sequences:
            literal_sequences_joined = ",".join(inputs.reference.literal_sequences)
            args.append(f"literal={literal_sequences_joined}")

        if inputs.operations.k_trim != "f":
            args.append(f"ktrim={inputs.operations.k_trim}")

        if inputs.operations.k_mask != "f" and inputs.operations.k_mask != "t":
            args.append(f"kmask={inputs.operations.k_mask}")

        if inputs.operations.quality_trim != "f":
            args.append(f"qtrim={inputs.operations.quality_trim}")

        if inputs.operations.ignore_bad_quality:
            args.append("ignorebadquality")

        if inputs.operations.max_length:
            args.append(f"masklength={inputs.operations.max_length}")

        if inputs.header_parsing.barcode_files:
            barcodes = ",".join(barcodes)
            args.append(f"barcodes={barcodes}")

        process_outputs = Parallel(n_jobs=n_jobs)(
            run_bbduk(input_reads=input_set, bbduk_inputs=args, paired_end=True)
            for input_set in input_reads
        )

        for output in process_outputs:
            if output[0]:
                self.error("BBDuk failed.", output[1])

        self.progress(0.7)

        statistics = []
        for stats in input_reads:
            with open(stats["stats"], "rb") as orig_file:
                with gzip.open(f"{stats['stats']}.gz", "wb") as zipped_file:
                    zipped_file.writelines(orig_file)
            statistics.append(f"{stats['stats']}.gz")

        output_path = Path("./fastqc")
        output_path.mkdir(exist_ok=True)

        rename_preprocessed_files(input_files=input_reads, paired_end=True)

        fastqgz = [fastq["original_name"] for fastq in input_reads]
        fastqgz2 = [fastq["original_name2"] for fastq in input_reads]
        fastqc_inputs = fastqgz + ["--extract", f"--outdir={str(output_path)}"]
        fastqc2_inputs = fastqgz2 + ["--extract", f"--outdir={str(output_path)}"]

        if inputs.fastqc.nogroup:
            fastqc_inputs.append("--no-group")
            fastqc2_inputs.append("--no-group")

        return_code, _, stderr = Cmd["fastqc"][fastqc_inputs] & TEE(retcode=None)
        if return_code:
            self.error("FastQC failed. ", stderr)
        return_code, _, stderr = Cmd["fastqc"][fastqc2_inputs] & TEE(retcode=None)
        if return_code:
            self.error("FastQC failed. ", stderr)

        for fastqc_zip in output_path.glob("*_fastqc.zip"):
            shutil.move(str(fastqc_zip), ".")

        mate1_fastqc, mate1_fastqc_url = prepare_fastqc(
            fastqgz=fastqgz,
            error=self.error,
        )

        mate2_fastqc, mate2_fastqc_url = prepare_fastqc(
            fastqgz=fastqgz2,
            error=self.error,
        )

        if (
            inputs.operations.k_trim == "f"
            and inputs.operations.quality_trim == "f"
            and not inputs.reference.sequences
            and not inputs.reference.literal_sequences
            and inputs.operations.min_average_quality <= 0
        ):
            self.warning(
                "Reference sequences, trimming mode, and minimum average quality are "
                "unspecified. Only filtering of reads by length is applied."
            )

        outputs.fastq = fastqgz
        outputs.fastq2 = fastqgz2
        outputs.fastqc_url = mate1_fastqc_url
        outputs.fastqc_url2 = mate2_fastqc_url
        outputs.fastqc_archive = mate1_fastqc
        outputs.fastqc_archive2 = mate2_fastqc
        outputs.statistics = statistics

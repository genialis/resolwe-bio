"""Run BBDuk."""
import gzip
import shutil
from pathlib import Path

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
    version = "2.6.0"
    category = "Trim"
    data_name = '{{ reads|sample_name|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.9.0"}
        },
        "resources": {
            "cores": 10,
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
        show_advanced = BooleanField(
            label="Show advanced parameters",
            default=False,
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

        reference = GroupField(Reference, label="Reference", hidden="!show_advanced")

        processing = GroupField(
            Processing, label="Processing parameters", hidden="!show_advanced"
        )

        operations = GroupField(
            Operations,
            label="Trimming, filtering and masking parameters.",
            hidden="!show_advanced",
        )

        header_parsing = GroupField(
            HeaderParsing, label="Header-parsing parameters", hidden="!show_advanced"
        )

        complexity = GroupField(
            Complexity, label="Complexity parameters", hidden="!show_advanced"
        )

        fastqc = GroupField(Fastqc, label="FastQC parameters", hidden="!show_advanced")

    class Output:
        """Output fields."""

        fastq = ListField(FileField(), label="Remaining reads")
        statistics = ListField(FileField(), label="Statistics")
        fastqc_url = ListField(FileHtmlField(), label="Quality control with FastQC")
        fastqc_archive = ListField(FileField(), label="Download FastQC archive")

    def run(self, inputs, outputs):
        """Run analysis."""
        reads_basename = Path(inputs.reads.output.fastq[0].path).name
        assert reads_basename.endswith(".fastq.gz")
        name = reads_basename[:-9]

        output_reads = "output_reads.fastq"
        output_statistics = "output_statistics.txt"
        input_reads = "input_reads.fastq.gz"
        input_references = "input_references.fasta.gz"
        input_barcodes = "input_barcodes.fasta.gz"
        final_reads = name + "_preprocessed.fastq"
        final_statistics = name + "_statistics.txt"

        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]]
            > input_reads
        )()
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
            f"in={input_reads}",
            f"out={output_reads}",
            f"stats={output_statistics}",
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
            f"threads={self.requirements.resources.cores}",
            f"-Xmx{int(0.85*self.requirements.resources.memory)}m",
        ]

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

        return_code, stdout, stderr = Cmd["bbduk.sh"][args] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error(f"Error occured with BBDuk:{stdout}\n{stderr}")

        self.progress(0.7)

        output_reads = Path(output_reads)
        output_statistics = Path(output_statistics)
        output_reads.rename(final_reads)
        output_statistics.rename(final_statistics)

        Cmd["pigz"][final_reads]()
        Cmd["pigz"][final_statistics]()

        self.progress(0.75)

        args_fastqc = [
            f"{final_reads}.gz",
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]

        if inputs.fastqc.nogroup:
            args_fastqc.append("--nogroup")

        return_code, stdout, stderr = Cmd["fastqc.sh"][args_fastqc] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error(f"Error while preparing FASTQC report:{stdout}\n{stderr}")

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

        outputs.fastq = [f"{final_reads}.gz"]
        outputs.statistics = [f"{final_statistics}.gz"]


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
    version = "2.6.0"
    category = "Trim"
    data_name = '{{ reads|sample_name|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.9.0"}
        },
        "resources": {
            "cores": 10,
            "memory": 8192,
        },
    }

    class Input:
        """Input fields to process BBDukSingle."""

        reads = DataField("reads:fastq:paired", label="Reads")
        min_length = IntegerField(
            label="Minimum length",
            default=10,
            description="Reads shorter than the minimum length will be discarded after trimming.",
        )
        show_advanced = BooleanField(
            label="Show advanced parameters",
            default=False,
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

        reference = GroupField(Reference, label="Reference", hidden="!show_advanced")

        processing = GroupField(
            Processing, label="Processing parameters", hidden="!show_advanced"
        )

        operations = GroupField(
            Operations,
            label="Trimming, filtering and masking parameters.",
            hidden="!show_advanced",
        )

        header_parsing = GroupField(
            HeaderParsing, label="Header-parsing parameters", hidden="!show_advanced"
        )

        complexity = GroupField(
            Complexity, label="Complexity parameters", hidden="!show_advanced"
        )

        fastqc = GroupField(Fastqc, label="FastQC parameters", hidden="!show_advanced")

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
        reads_mate1_basename = Path(inputs.reads.output.fastq[0].path).name
        assert reads_mate1_basename.endswith(".fastq.gz")
        name_mate1 = reads_mate1_basename[:-9]
        reads_mate2_basename = Path(inputs.reads.output.fastq2[0].path).name
        assert reads_mate2_basename.endswith(".fastq.gz")
        name_mate2 = reads_mate2_basename[:-9]

        input_mate1_reads = "input_mate1_reads.fastq.gz"
        input_mate2_reads = "input_mate2_reads.fastq.gz"
        input_references = "input_references.fasta.gz"
        input_barcodes = "input_barcodes.fasta.gz"
        output_mate1_reads = "output_mate1_reads.fastq"
        output_mate2_reads = "output_mate2_reads.fastq"
        output_statistics = "output_statistics.txt"
        final_mate1_reads = name_mate1 + "_preprocessed.fastq"
        final_mate2_reads = name_mate2 + "_preprocessed.fastq"
        final_statistics = name_mate1 + "_statistics.txt"

        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq]]
            > input_mate1_reads
        )()
        (
            Cmd["cat"][[reads.path for reads in inputs.reads.output.fastq2]]
            > input_mate2_reads
        )()
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
            f"in={input_mate1_reads}",
            f"in2={input_mate2_reads}",
            f"out={output_mate1_reads}",
            f"out2={output_mate2_reads}",
            f"stats={output_statistics}",
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
            f"threads={self.requirements.resources.cores}",
            f"-Xmx{int(0.85*self.requirements.resources.memory)}m",
        ]

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

        return_code, stdout, stderr = Cmd["bbduk.sh"][args] & TEE(retcode=None)
        if return_code:
            print(stdout, stderr)
            self.error(f"Error occured with BBDuk:{stdout}\n{stderr}")

        self.progress(0.7)

        output_mate1_reads = Path(output_mate1_reads)
        output_mate2_reads = Path(output_mate2_reads)
        output_statistics = Path(output_statistics)
        output_mate1_reads.rename(final_mate1_reads)
        output_mate2_reads.rename(final_mate2_reads)
        output_statistics.rename(final_statistics)

        Cmd["pigz"][final_mate1_reads]()
        Cmd["pigz"][final_mate2_reads]()
        Cmd["pigz"][final_statistics]()

        self.progress(0.75)

        args_fastqc1 = [
            f"{final_mate1_reads}.gz",
            "fastqc",
            "fastqc_archive",
            "fastqc_url",
        ]

        args_fastqc2 = [
            f"{final_mate2_reads}.gz",
            "fastqc",
            "fastqc_archive2",
            "fastqc_url2",
        ]

        for args in [args_fastqc1, args_fastqc2]:
            if inputs.fastqc.nogroup:
                args.append("--nogroup")

            return_code, stdout, stderr = Cmd["fastqc.sh"][args] & TEE(retcode=None)
            if return_code:
                print(stdout, stderr)
                self.error(f"Error while preparing FASTQC report:{stdout}\n{stderr}")

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

        outputs.fastq = [f"{final_mate1_reads}.gz"]
        outputs.fastq2 = [f"{final_mate2_reads}.gz"]
        outputs.statistics = [f"{final_statistics}.gz"]

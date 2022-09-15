"""RNA-seq Variant Calling Workflow (Beta)."""
from resolwe.process import (
    BooleanField,
    Data,
    DataField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    StringField,
)
from resolwe.process.models import Process as BioProcess


class WorkflowRnaseqVariantCalling(Process):
    """(Beta version) Identify variants in RNA-seq data.

    This pipeline follows GATK best practices recommendantions for
    variant calling with RNA-seq data.

    The pipeline steps include read alignment (STAR), data cleanup
    (MarkDuplicates), splitting reads that contain Ns in their cigar string
    (SplitNCigarReads), base quality recalibration (BaseRecalibrator,
    ApplyBQSR), variant calling (HaplotypeCaller) and variant filtering
    (VariantFiltration, SelectVariants).
    """

    slug = "workflow-rnaseq-variantcalling-beta"
    name = "RNA-seq Variant Calling (beta)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/dnaseq:6.3.1",
            },
        },
    }
    data_name = "RNA-seq Variants ({{ reads|name|default('?') }})"
    version = "1.0.0"
    process_type = "data:workflow:rnaseq:variants"
    category = "Pipeline"
    entity = {
        "type": "sample",
    }

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Input sample (FASTQ)",
            description="Input data in FASTQ format.",
        )
        ref_seq = DataField(
            data_type="seq:nucleotide", label="Reference FASTA sequence"
        )
        genome = DataField(
            data_type="index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
        )
        preprocessing = BooleanField(
            label="Perform reads processing with BBDuk",
            default=True,
            description="If your reads have not been processed, set this to True.",
        )
        dbsnp = DataField(
            data_type="variants:vcf",
            label="dbSNP file",
            description="File with known variants.",
        )
        indels = ListField(
            DataField(data_type="variants:vcf"),
            label="Known INDEL sites",
            required=False,
        )
        intervals = DataField(
            data_type="bed",
            label="Intervals (from BED file)",
            description="Use this option to perform the analysis over only part "
            "of the genome.",
            required=False,
        )
        read_group = StringField(
            label="Replace read groups in BAM",
            description="Replace read groups in a BAM file. This argument enables the "
            "user to replace all read groups in the INPUT file with a single new read "
            "group and assign all reads to this read group in the OUTPUT BAM file. "
            "Addition or replacement is performed using Picard's AddOrReplaceReadGroups "
            "tool. Input should take the form of -name=value delimited by a "
            '";", e.g. "-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1". '
            "See tool's documentation for more information on tag names. Note that "
            "PL, LB, PU and SM are required fields. See caveats of rewriting read groups "
            "in the documentation.",
            default="-ID=1;-LB=GENIALIS;-PL=ILLUMINA;-PU=BARCODE;-SM=SAMPLENAME1",
        )
        filter_expressions = ListField(
            StringField(),
            label="Expressions used with INFO fields to filter",
            description="VariantFiltration accepts any number of JEXL expressions "
            "(so you can have two named filters by using --filter-name One "
            "--filter-expression 'X < 1' --filter-name Two --filter-expression 'X > 2'). "
            "It is preferable to use multiple expressions, each specifying an individual "
            "filter criteria, to a single compound expression that specifies multiple "
            "filter criteria. Input expressions one by one and press ENTER after each "
            "expression. Examples of filter expression: 'FS > 30', 'DP > 10'.",
            default=["FS > 30.0", "QD < 2.0", "DP < 10.0"],
        )
        filter_name = ListField(
            StringField(),
            label="Names to use for the list of filters",
            description="This name is put in the FILTER field for variants that get "
            "filtered. Note that there must be a 1-to-1 mapping between filter expressions "
            "and filter names. Input expressions one by one and press ENTER after each name. "
            "Warning: filter names should be in the same order as filter expressions. "
            "Example: you specified filter expressions 'FS > 30' and 'DP > 10', now "
            "specify filter names 'FS' and 'DP'.",
            default=["FS", "QD", "DP"],
        )
        exclude_filtered = BooleanField(
            label="Exclude filtered sites",
            default=True,
            description="If this flag is enabled, sites that have been marked as filtered "
            "(i.e. have anything other than `.` or `PASS` in the FILTER field) will be excluded "
            "from the output.",
        )

        class Bbduk:
            """Preprocessing options (BBDuk)."""

            adapters = ListField(
                DataField(data_type="seq:nucleotide"),
                label="Adapters",
                required=False,
                description="Provide a list of sequencing adapters files (.fasta) "
                "to be removed by BBDuk.",
            )
            custom_adapter_sequences = ListField(
                StringField(),
                label="Custom adapter sequences",
                required=False,
                default=[],
                description="Custom adapter sequences can be specified by inputting them "
                "one by one and pressing Enter after each sequence.",
            )
            kmer_length = IntegerField(
                label="K-mer length [k=]",
                default=23,
                description="Kmer length used for finding contaminants. "
                "Contaminants shorter than kmer length will not be found. "
                "Kmer length must be at least 1.",
            )
            min_k = IntegerField(
                label="Minimum k-mer length at right end of reads used for trimming [mink=]",
                default=11,
                disabled="bbduk.adapters.length === 0 && bbduk.custom_adapter_sequences.length === 0",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for k-mers [hammingdistance=]",
                default=1,
                description="Hamming distance i.e. the number of mismatches allowed in the kmer.",
            )
            maxns = IntegerField(
                label="Max Ns after trimming [maxns=]",
                default=-1,
                description="If non-negative, reads with more Ns than this (after trimming) will be discarded.",
            )
            trim_quality = IntegerField(
                label="Average quality below which to trim region [trimq=]",
                default=28,
                description="Phred algorithm is used, which is more accurate than naive trimming.",
            )
            min_length = IntegerField(
                label="Minimum read length [minlength=]",
                default=30,
                description="Reads shorter than minimum read length after trimming are discarded.",
            )
            quality_encoding_offset = StringField(
                label="Quality encoding offset",
                choices=[
                    ("33", "Sanger / Illumina 1.8+"),
                    ("64", "Illumina up to 1.3+, 1.5+"),
                    ("auto", "Auto"),
                ],
                default="auto",
                description="Quality encoding offset for input FASTQ files.",
            )
            ignore_bad_quality = BooleanField(
                label="Ignore bad quality",
                default=False,
                description="Don't crash if quality values appear to be incorrect.",
            )

        class Alignment:
            """Alignment options (STAR)."""

            two_pass_mode = BooleanField(
                label="Use two pass mode [--twopassMode]",
                default=True,
                description="Use two-pass maping instead of first-pass only. In two-pass mode we "
                "first perform first-pass mapping, extract junctions, insert them into genome "
                "index, and re-map all reads in the second mapping pass.",
            )
            out_unmapped = BooleanField(
                label="Output unmapped reads (SAM) [--outSAMunmapped Within]",
                default=True,
                description="Output of unmapped reads in the SAM format.",
            )

        class HaplotypeCaller:
            """GATK HaplotypeCaller options."""

            stand_call_conf = IntegerField(
                label="Min call confidence threshold",
                default=20,
                description="The minimum phred-scaled confidence threshold at which "
                "variants should be called.",
            )
            soft_clipped = BooleanField(
                label="Do not analyze soft clipped bases in the reads",
                default=True,
                description="Suitable option for RNA-seq variant calling.",
            )

        class SelectVariants:
            """GATK SelectVariants options."""

            select_type = ListField(
                StringField(),
                label="Select only a certain type of variants from the input file",
                description="This argument selects particular kinds of variants out of a list. If "
                "left empty, there is no type selection and all variant types are considered for "
                "other selection criteria. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, "
                "NO_VARIATION. Can be specified multiple times.",
                default=["SNP"],
            )

        class Advanced:
            """Advanced options."""

            java_gc_threads = IntegerField(
                label="Java ParallelGCThreads",
                default=2,
                description="Sets the number of threads used during parallel phases of "
                "the garbage collectors.",
            )
            max_heap_size = IntegerField(
                label="Java maximum heap size (Xmx)",
                default=12,
                description="Set the maximum Java heap size (in GB).",
            )

        bbduk = GroupField(
            Bbduk, label="Preprocessing with BBDuk", hidden="!preprocessing"
        )
        alignment = GroupField(
            Alignment,
            label="Alignment with STAR",
        )
        haplotype_caller = GroupField(
            HaplotypeCaller, label="Options for HaplotypeCaller"
        )
        select_variants = GroupField(
            SelectVariants,
            label="Options for GATK SelectVariants",
            hidden="!exclude_filtered",
        )
        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields."""

        # Workflow does not have outputs.

    def run(self, inputs, outputs):
        """Run the workflow."""

        if (
            inputs.ref_seq.output.species != inputs.genome.output.species
            or inputs.dbsnp.output.species != inputs.genome.output.species
        ):
            self.error(
                "All input files must be from the same species. "
                f"FASTA reference sequence is from {inputs.ref_seq.output.species}, "
                f"STAR index is from {inputs.genome.output.species} and "
                f"DBSNP file is from {inputs.dbsnp.output.species}"
            )

        if (
            inputs.ref_seq.output.build != inputs.genome.output.build
            or inputs.dbsnp.output.build != inputs.genome.output.build
        ):
            self.error(
                "All input files must have the same build. "
                f"FASTA reference sequence is based on {inputs.ref_seq.output.build}, "
                f"STAR index is based on {inputs.genome.output.build} and "
                f"DBSNP file is based on {inputs.dbsnp.output.build}. "
            )

        if inputs.intervals:
            if inputs.intervals.output.species != inputs.genome.output.species:
                self.error(
                    "All input files must be from the same species. "
                    f"Intervals BED file is from {inputs.intervals.output.species}, "
                    f"while STAR index is from {inputs.genome.output.species}."
                )
            if inputs.intervals.output.build != inputs.genome.output.build:
                self.error(
                    "All input files must have the same build. "
                    f"Intervals BED file is based on {inputs.intervals.output.build}, "
                    f"while STAR index is based on {inputs.genome.output.build}."
                )

        if inputs.indels:
            for indel in inputs.indels:
                if indel.output.species != inputs.genome.output.species:
                    self.error(
                        "All input files must be from the same species. "
                        f"File with INDELs is from {indel.output.species}, "
                        f"while STAR index is from {inputs.genome.output.species}."
                    )
                if indel.output.build != inputs.genome.output.build:
                    self.error(
                        "All input file must have the same build. "
                        f"File with INDELs is based on {indel.output.build}, "
                        f"while STAR index is based on {inputs.genome.output.build}."
                    )

        if inputs.preprocessing:

            # BBDuk options used are the same as the default options
            # in workflow-bbduk-star-featurecounts-qc
            input_preprocessing = {
                "reads": inputs.reads,
                "min_length": inputs.bbduk.min_length,
                "reference": {
                    "sequences": inputs.bbduk.adapters or [],
                    "literal_sequences": inputs.bbduk.custom_adapter_sequences,
                },
                "processing": {
                    "kmer_length": inputs.bbduk.kmer_length,
                    "hamming_distance": inputs.bbduk.hamming_distance,
                },
                "operations": {
                    "quality_trim": "r",
                    "trim_quality": inputs.bbduk.trim_quality,
                    "quality_encoding_offset": inputs.bbduk.quality_encoding_offset,
                    "ignore_bad_quality": inputs.bbduk.ignore_bad_quality,
                    "maxns": 1,
                },
            }

            if inputs.bbduk.adapters or inputs.bbduk.custom_adapter_sequences:
                input_preprocessing["operations"]["k_trim"] = "r"
            else:
                input_preprocessing["operations"]["k_trim"] = "f"

            if inputs.bbduk.adapters or inputs.bbduk.custom_adapter_sequences:
                input_preprocessing["operations"]["min_k"] = inputs.bbduk.min_k
            else:
                input_preprocessing["operations"]["min_k"] = -1

            if inputs.reads.type.startswith("data:reads:fastq:single:"):
                slug_bbduk = "bbduk-single"
            elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
                slug_bbduk = "bbduk-paired"
            else:
                self.error("Wrong reads input type.")

            preprocessing = Data.create(
                process=BioProcess.get_latest(slug=slug_bbduk),
                input=input_preprocessing,
                name=f"Trimmed ({inputs.reads.name})",
            )

        input_star = {
            "reads": preprocessing if inputs.preprocessing else inputs.reads,
            "genome": inputs.genome,
            "two_pass_mapping": {
                "two_pass_mode": inputs.alignment.two_pass_mode,
            },
            "output_options": {
                "out_unmapped": inputs.alignment.out_unmapped,
            },
        }

        alignment = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input=input_star,
            name=f"Aligned ({inputs.reads.name})",
        )

        input_preprocess = {
            "bam": alignment,
            "ref_seq": inputs.ref_seq,
            "known_sites": [inputs.dbsnp],
            "read_group": inputs.read_group,
        }
        if inputs.indels:
            for indel in inputs.indels:
                input_preprocess["known_sites"].append(indel)

        preprocess = Data.create(
            process=BioProcess.get_latest(slug="rnaseq-vc-preprocess"),
            input=input_preprocess,
            name=f"Preprocessed ({inputs.reads.name})",
        )

        input_hc = {
            "alignment": preprocess,
            "genome": inputs.ref_seq,
            "dbsnp": inputs.dbsnp,
            "stand_call_conf": inputs.haplotype_caller.stand_call_conf,
            "mbq": 10,
            "advanced": {
                "soft_clipped": inputs.haplotype_caller.soft_clipped,
                "java_gc_threads": inputs.advanced.java_gc_threads,
                "max_heap_size": inputs.advanced.max_heap_size,
            },
        }

        if inputs.intervals:
            input_hc["intervals_bed"] = inputs.intervals

        hc = Data.create(
            process=BioProcess.get_latest(slug="vc-gatk4-hc"),
            input=input_hc,
            name=f"Genotypes ({inputs.reads.name})",
        )

        input_filtration = {
            "vcf": hc,
            "ref_seq": inputs.ref_seq,
            "filter_expressions": inputs.filter_expressions,
            "filter_name": inputs.filter_name,
            "advanced": {
                "cluster": 3,
                "window": 35,
                "java_gc_threads": inputs.advanced.java_gc_threads,
                "max_heap_size": inputs.advanced.max_heap_size,
            },
        }

        filtration = Data.create(
            process=BioProcess.get_latest(slug="gatk-variant-filtration-single"),
            input=input_filtration,
            name=f"Filtered genotypes ({inputs.reads.name})",
        )

        if inputs.exclude_filtered:
            input_select = {
                "vcf": filtration,
                "select_type": inputs.select_variants.select_type,
                "exclude_filtered": True,
                "advanced_options": {
                    "java_gc_threads": inputs.advanced.java_gc_threads,
                    "max_heap_size": inputs.advanced.max_heap_size,
                },
            }

            Data.create(
                process=BioProcess.get_latest(slug="gatk-select-variants-single"),
                input=input_select,
                name=f"Selected genotypes ({inputs.reads.name})",
            )

        multiqc_tools = [inputs.reads, alignment]

        if inputs.preprocessing:
            multiqc_tools.append(preprocessing)

        input_multiqc = {"data": multiqc_tools}

        Data.create(
            process=BioProcess.get_latest(slug="multiqc"),
            input=input_multiqc,
        )

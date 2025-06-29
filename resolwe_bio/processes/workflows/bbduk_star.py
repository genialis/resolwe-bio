"""STAR-based RNA-seq pipeline."""

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


class WorkflowSTAR(Process):
    """STAR-based RNA-seq pipeline.

    First, reads are preprocessed by __BBDuk__ which removes adapters, trims
    reads for quality from the 3'-end, and discards reads that are too short
    after trimming. Compared to similar tools, BBDuk is regarded for its
    computational efficiency. Next, preprocessed reads are aligned by __STAR__
    aligner. At the time of implementation, STAR is considered a
    state-of-the-art tool that consistently produces accurate results from
    diverse sets of reads, and performs well even with default settings. STAR aligner
    counts and reports the number of aligned reads per gene while mapping.
    STAR version used is 2.7.10b. For more information see
    [this comparison of RNA-seq aligners](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/).

    rRNA contamination rate in the sample is determined using the STAR aligner.
    Quality-trimmed reads are downsampled (using __Seqtk__ tool) and aligned to the
    rRNA reference sequences. The alignment rate indicates the percentage of the
    reads in the sample that are derived from the rRNA sequences. Final step of the
    workflow is QoRTs QC analysis with downsampled reads.
    """

    slug = "workflow-bbduk-star-qc"
    name = "STAR-based gene quantification workflow"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    version = "1.8.0"
    entity = {
        "type": "sample",
    }
    process_type = "data:workflow:rnaseq:star:qc"
    category = "Pipeline"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Reads (FASTQ)",
            description="Reads in FASTQ file, single or paired end.",
        )
        genome = DataField(
            data_type="index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
        )
        annotation = DataField(
            data_type="annotation",
            label="Annotation",
            description="GTF and GFF3 annotation formats are supported.",
        )
        assay_type = StringField(
            label="Assay type",
            choices=[
                ("non_specific", "Strand non-specific"),
                ("forward", "Strand-specific forward"),
                ("reverse", "Strand-specific reverse"),
                ("auto", "Detect automatically"),
            ],
            description="In strand non-specific assay a read is considered overlapping with a "
            "feature regardless of whether it is mapped to the same or the opposite "
            "strand as the feature. In strand-specific forward assay and single "
            "reads, the read has to be mapped to the same strand as the feature. "
            "For paired-end reads, the first read has to be on the same strand and "
            "the second read on the opposite strand. In strand-specific reverse "
            "assay these rules are reversed.",
            default="non_specific",
        )
        cdna_index = DataField(
            data_type="index:salmon",
            label="Indexed cDNA reference sequence",
            required=False,
            description="Transcriptome index file created using the Salmon indexing tool. "
            "cDNA (transcriptome) sequences used for index file creation must be "
            "derived from the same species as the input sequencing reads to "
            "obtain the reliable analysis results.",
            hidden="assay_type != 'auto'",
        )
        rrna_reference = DataField(
            data_type="index:star",
            label="Indexed rRNA reference sequence",
            description="Reference sequence index prepared by STAR aligner indexing tool.",
        )
        globin_reference = DataField(
            data_type="index:star",
            label="Indexed Globin reference sequence",
            description="Reference sequence index prepared by STAR aligner indexing tool.",
        )
        qc_evaluation = BooleanField(
            label="Enable QC evaluation in MultiQC",
            default=True,
            description="Control whether QC evaluation and annotation is performed in the MultiQC process.",
        )

        class Preprocessing:
            """Preprocessing with BBDuk."""

            adapters = ListField(
                inner=DataField(data_type="seq:nucleotide"),
                label="Adapters",
                required=False,
                description="FASTA file(s) with adapters.",
            )
            custom_adapter_sequences = ListField(
                inner=StringField(),
                label="Custom adapter sequences",
                required=False,
                default=[],
                description="Custom adapter sequences can be specified by inputting them "
                "one by one and pressing Enter after each sequence.",
            )
            kmer_length = IntegerField(
                label="K-mer length [k=]",
                default=23,
                description="K-mer length used for finding contaminants. "
                "Contaminants shorter than k-mer length will not be found. "
                "K-mer length must be at least 1.",
            )
            min_k = IntegerField(
                label="Minimum k-mer length at right end of reads used for trimming [mink=]",
                default=11,
                disabled="preprocessing.adapters.length === 0 && preprocessing.custom_adapter_sequences.length === 0",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for k-mers [hammingdistance=]",
                default=1,
                description="Hamming distance i.e. the number of mismatches allowed in the k-mer.",
            )
            maxns = IntegerField(
                label="Max Ns after trimming [maxns=]",
                default=-1,
                description="If non-negative, reads with more Ns than this (after trimming) will be discarded.",
            )
            trim_quality = IntegerField(
                label="Average quality below which to trim region [trimq=]",
                default=10,
                description="Phred algorithm is used, which is more accurate than naive trimming.",
            )
            min_length = IntegerField(
                label="Minimum read length [minlength=]",
                default=20,
                description="Reads shorter than minimum read length after trimming are discarded.",
            )
            quality_encoding_offset = StringField(
                label="Quality encoding offset [qin=]",
                choices=[
                    ("33", "Sanger / Illumina 1.8+"),
                    ("64", "Illumina up to 1.3+, 1.5+"),
                    ("auto", "Auto"),
                ],
                default="auto",
                description="Quality encoding offset for input FASTQ files.",
            )
            ignore_bad_quality = BooleanField(
                label="Ignore bad quality [ignorebadquality]",
                default=False,
                description="Don't crash if quality values appear to be incorrect.",
            )

        class Alignment:
            """Alignment with STAR."""

            unstranded = BooleanField(
                label="The data is unstranded [--outSAMstrandField intronMotif]",
                default=False,
                description="For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced "
                "alignments with XS strand attribute, which STAR will generate with "
                "--outSAMstrandField intronMotif option. As required, the XS strand "
                "attribute will be generated for all alignments that contain splice "
                "junctions. The spliced alignments that have undefined strand "
                "(i.e. containing only non-canonical unannotated junctions) will be "
                "suppressed. If you have stranded RNA-seq data, you do not need to "
                "use any specific STAR options. Instead, you need to run Cufflinks with "
                "the library option --library-type options. For example, "
                "cufflinks --library-type fr-firststrand should be used for the standard "
                "dUTP protocol, including Illumina's stranded Tru-Seq. "
                "This option has to be used only for Cufflinks runs and not for STAR runs.",
            )
            noncannonical = BooleanField(
                label="Remove non-canonical junctions (Cufflinks compatibility)",
                default=False,
                description="It is recommended to remove the non-canonical junctions for Cufflinks "
                "runs using --outFilterIntronMotifs RemoveNoncanonical.",
            )

            class ChimericReadsOptions:
                """Chimeric reads options."""

                chimeric = BooleanField(
                    label="Detect chimeric and circular alignments [--chimOutType]",
                    default=False,
                    description="Detect chimeric and circular alignments. If set to True, "
                    "certain STAR parameters are set to values recommended for maximum sensitivity "
                    "by Arriba documentation.",
                )

                chim_out_type = StringField(
                    label="Chimeric output type [--chimOutType]",
                    default="WithinBAM HardClip",
                    choices=[
                        ("WithinBAM", "WithinBAM"),
                        ("WithinBAM HardClip", "WithinBAM HardClip"),
                        ("WithinBAM SoftClip", "WithinBAM SoftClip"),
                        ("SeparateSAMold", "SeparateSAMold"),
                    ],
                    description="Type of chimeric output produced by STAR. If WithinBAM is "
                    "selected, 'WithinBAM HardClip' in the CIGAR is used.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_segment_min = IntegerField(
                    label="Minimum length of chimeric segment [--chimSegmentMin]",
                    default=10,
                    description="Minimum length of chimeric segment to be reported. "
                    "Each chimeric alignment consists of two segments. Each segment is non-chimeric on "
                    "its own, but the segments are chimeric to each other (i.e. the segments belong "
                    "to different chromosomes, or different strands, or are far from each other). "
                    "Both segments may contain splice junctions, and one of the segments may contain "
                    "portions of both mates. --chimSegmentMin parameter controls the minimum mapped "
                    "length of the two segments that is allowed. For example, if you have 2x75 reads "
                    "and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome "
                    "and 20b on the other will be output, while 135 + 15 won't be. Additionally, "
                    "filtering_options.out_multimap_max is set to 50 if detect_chimeric is set to True.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_junction_overhang_min = IntegerField(
                    label="Minimum overhang for chimeric junctions [--chimJunctionOverhangMin]",
                    default=10,
                    description="Minimum overhang for chimeric junctions.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_score_drop_max = IntegerField(
                    label="Maximum drop in chimeric score [--chimScoreDropMax]",
                    default=30,
                    description="Maximum drop in chimeric score for chimeric alignments "
                    "relative to the read length.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_score_junction_non_gtag = IntegerField(
                    label="Score for non-GTAG chimeric junctions [--chimScoreJunctionNonGTAG]",
                    default=0,
                    description="Penalty for chimeric junctions if they are non-GTAG.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_score_separation = IntegerField(
                    label="Minimum score separation between chimeric alignments [--chimScoreSeparation]",
                    default=1,
                    description="Minimum score separation between the best chimeric alignment "
                    "and the next best.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_segment_read_gap_max = IntegerField(
                    label="Maximum read gap for chimeric segments [--chimSegmentReadGapMax]",
                    default=3,
                    description="Maximum gap in the read sequence between chimeric segments.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                chim_multimap_nmax = IntegerField(
                    label="Maximum number of chimeric multi-mapping segments [--chimMultimapNmax]",
                    default=50,
                    description="Maximum number of chimeric alignments reported for multimapped reads.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                pe_overlap_n_bases_min = IntegerField(
                    label="Minimum number of overlapping bases in PE [--peOverlapNbasesMin]",
                    default=10,
                    description="Minimum number of overlapping bases to trigger merging mate pair.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                align_spliced_mate_map_lmin_over_lmate = FloatField(
                    label="Minimum mapped length of spliced mates over mate "
                    "length [--alignSplicedMateMapLminOverLmate]",
                    default=0.5,
                    description="Minimum mapped length of the spliced mate, divided by the mate length.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

                align_sj_stitch_mismatch_n_max = ListField(
                    IntegerField(),
                    label="Maximum number of mismatches for stitching SJ [--alignSJstitchMismatchNmax]",
                    default=[5, -1, 5, 5],
                    description="Maximum number of mismatches that are allowed in the seed region for "
                    "stitching splice junctions. Enter four integers.",
                    disabled="!alignment.chimeric_reads.chimeric",
                )

            class TranscriptOutputOptions:
                """Transcript coordinate output options."""

                quant_mode = BooleanField(
                    label="Output in transcript coordinates [--quantMode]",
                    default=False,
                    description="With --quantMode TranscriptomeSAM option STAR will output alignments "
                    "translated into transcript coordinates in the Aligned.toTranscriptome.out.bam "
                    "file (in addition to alignments in genomic coordinates in Aligned.*.sam/bam "
                    "files). These transcriptomic alignments can be used with various transcript "
                    "quantification software that require reads to be mapped to transcriptome, such "
                    "as RSEM or eXpress.",
                )
                single_end = BooleanField(
                    label="Allow soft-clipping and indels [--quantTranscriptomeBan Singleend]",
                    default=False,
                    disabled="!t_coordinates.quant_mode",
                    description="By default, the output satisfies RSEM requirements: soft-clipping or "
                    "indels are not allowed. Use --quantTranscriptomeBan Singleend to allow "
                    "insertions, deletions and soft-clips in the transcriptomic alignments, which "
                    "can be used by some expression quantification softwares (e.g. eXpress).",
                )

            class FilteringOptions:
                """Output filtering options."""

                out_filter_type = StringField(
                    label="Type of filtering [--outFilterType]",
                    default="Normal",
                    choices=[
                        ("Normal", "Normal"),
                        ("BySJout", "BySJout"),
                    ],
                    description="Normal: standard filtering using only current alignment; BySJout: "
                    "keep only those reads that contain junctions that passed filtering into "
                    "SJ.out.tab.",
                )
                out_multimap_max = IntegerField(
                    label="Maximum number of loci [--outFilterMultimapNmax]",
                    required=False,
                    description="Maximum number of loci the read is allowed to map to. Alignments "
                    "(all of them) will be output only if the read maps to no more loci than this "
                    "value. Otherwise no alignments will be output, and the read will be counted as "
                    "'mapped to too many loci' (default: 10).",
                )
                out_mismatch_max = IntegerField(
                    label="Maximum number of mismatches [--outFilterMismatchNmax]",
                    required=False,
                    description="Alignment will be output only if it has fewer mismatches than this "
                    "value (default: 10). Large number (e.g. 999) switches off this filter.",
                )
                out_mismatch_nl_max = FloatField(
                    label="Maximum no. of mismatches (map length) [--outFilterMismatchNoverLmax]",
                    required=False,
                    range=[0.0, 1.0],
                    description="Alignment will be output only if its ratio of mismatches to *mapped* "
                    "length is less than or equal to this value (default: 0.3). The value should be "
                    "between 0.0 and 1.0.",
                )
                out_score_min = IntegerField(
                    label="Minimum alignment score [--outFilterScoreMin]",
                    required=False,
                    description="Alignment will be output only if its score is higher than or equal "
                    "to this value (default: 0).",
                )
                out_mismatch_nrl_max = FloatField(
                    label="Maximum no. of mismatches (read length) [--outFilterMismatchNoverReadLmax]",
                    required=False,
                    range=[0.0, 1.0],
                    description="Alignment will be output only if its ratio of mismatches to *read* "
                    "length is less than or equal to this value (default: 1.0). Using 0.04 for "
                    "2x100bp, the max number of mismatches is calculated as 0.04*200=8 for the paired "
                    "read. The value should be between 0.0 and 1.0.",
                )

            class AlignmentOptions:
                """Alignment and Seeding."""

                align_overhang_min = IntegerField(
                    label="Minimum overhang [--alignSJoverhangMin]",
                    required=False,
                    description="Minimum overhang (i.e. block size) for spliced alignments "
                    "(default: 5).",
                )
                align_sjdb_overhang_min = IntegerField(
                    label="Minimum overhang (sjdb) [--alignSJDBoverhangMin]",
                    required=False,
                    description="Minimum overhang (i.e. block size) for annotated (sjdb) spliced "
                    "alignments (default: 3).",
                )
                align_intron_size_min = IntegerField(
                    label="Minimum intron size [--alignIntronMin]",
                    required=False,
                    description="Minimum intron size: the genomic gap is considered an intron if its "
                    "length >= alignIntronMin, otherwise it is considered Deletion (default: 21).",
                )
                align_intron_size_max = IntegerField(
                    label="Maximum intron size [--alignIntronMax]",
                    required=False,
                    description="Maximum intron size, if 0, max intron size will be determined by "
                    "(2pow(winBinNbits)*winAnchorDistNbins)(default: 0).",
                )
                align_gap_max = IntegerField(
                    label="Minimum gap between mates [--alignMatesGapMax]",
                    required=False,
                    description="Maximum gap between two mates, if 0, max intron gap will be "
                    "determined by (2pow(winBinNbits)*winAnchorDistNbins) (default: 0).",
                )
                align_end_alignment = StringField(
                    label="Read ends alignment [--alignEndsType]",
                    choices=[
                        ("Local", "Local"),
                        ("EndToEnd", "EndToEnd"),
                        ("Extend5pOfRead1", "Extend5pOfRead1"),
                        ("Extend5pOfReads12", "Extend5pOfReads12"),
                    ],
                    description="Type of read ends alignment (default: Local). Local: standard local "
                    "alignment with soft-clipping allowed. EndToEnd: force end-to-end read alignment, "
                    "do not soft-clip. Extend5pOfRead1: fully extend only the 5p of the read1, all "
                    "other ends: local alignment. Extend5pOfReads12: fully extend only the 5' of the "
                    "both read1 and read2, all other ends use local alignment.",
                    default="Local",
                )

            class TwoPassOptions:
                """Two-pass mapping options."""

                two_pass_mode = BooleanField(
                    label="Use two pass mode [--twopassMode]",
                    default=True,
                    description="Use two-pass maping instead of first-pass only. In two-pass mode we "
                    "first perform first-pass mapping, extract junctions, insert them into genome "
                    "index, and re-map all reads in the second mapping pass.",
                )

            class OutputOptions:
                """Output options."""

                out_unmapped = BooleanField(
                    label="Output unmapped reads (SAM) [--outSAMunmapped Within]",
                    default=True,
                    description="Output of unmapped reads in the SAM format.",
                )
                out_sam_attributes = StringField(
                    label="Desired SAM attributes [--outSAMattributes]",
                    default="Standard",
                    choices=[
                        ("Standard", "Standard"),
                        ("All", "All"),
                        ("NH HI NM MD", "NH HI NM MD"),
                        ("None", "None"),
                    ],
                    description="A string of desired SAM attributes, in the order desired for the "
                    "output SAM.",
                )
                out_rg_line = StringField(
                    label="SAM/BAM read group line [--outSAMattrRGline]",
                    required=False,
                    description="The first word contains the read group identifier and must start "
                    "with ID:, e.g. --outSAMattrRGline ID:xxx CN:yy ”DS:z z z” xxx will be added as "
                    "RG tag to each output alignment. Any spaces in the tag values have to be double "
                    "quoted. Comma separated RG lines corresponds to different (comma separated) input "
                    "files in -readFilesIn. Commas have to be surrounded by spaces, e.g. "
                    "-outSAMattrRGline ID:xxx , ID:zzz ”DS:z z” , ID:yyy DS:yyyy.",
                )

            chimeric_reads = GroupField(
                ChimericReadsOptions,
                label="Chimeric reads options",
            )
            transcript_output = GroupField(
                TranscriptOutputOptions,
                label="Transcript coordinate output options",
            )
            filtering_options = GroupField(
                FilteringOptions,
                label="Output filtering options",
            )
            alignment_options = GroupField(
                AlignmentOptions,
                label="Alignment options",
            )
            two_pass_mapping = GroupField(TwoPassOptions, label="Two-pass mapping")
            output_options = GroupField(
                OutputOptions,
                label="Output options",
            )

        class Quantification:
            """Quantification."""

            n_reads = IntegerField(
                label="Number of reads in subsampled alignment file for strandedness detection",
                default=5000000,
                hidden="assay_type != 'auto'",
                description="Alignment (.bam) file subsample size. Increase the number of reads "
                "to make automatic detection more reliable. Decrease the number of "
                "reads to make automatic detection run faster.",
            )

        class Downsampling:
            """Downsampling (Seqtk)."""

            n_reads = IntegerField(
                label="Number of reads",
                default=1000000,
                description="Number of reads to include in downsampling.",
            )

            class Advanced:
                """Advanced options for downsampling."""

                seed = IntegerField(
                    label="Seed [-s]",
                    default=11,
                    description="Using the same random seed makes reads downsampling more reproducible "
                    "in different environments.",
                )
                fraction = FloatField(
                    label="Fraction of reads used",
                    required=False,
                    range=[0.0, 1.0],
                    description="Use the fraction of reads [0.0 - 1.0] from the original input file instead "
                    "of the absolute number of reads. If set, this will override the 'Number of reads' "
                    "input parameter.",
                )
                two_pass = BooleanField(
                    label="2-pass mode [-2]",
                    default=False,
                    description="Enable two-pass mode when downsampling. Two-pass mode is twice "
                    "as slow but with much reduced memory.",
                )

            advanced = GroupField(
                Advanced,
                label="Advanced options for downsampling",
            )

        preprocessing = GroupField(
            Preprocessing,
            label="Preprocessing with BBDuk",
        )
        alignment = GroupField(
            Alignment,
            label="Alignment with STAR",
        )
        quantification = GroupField(
            Quantification,
            label="Quantification",
        )
        downsampling = GroupField(
            Downsampling,
            label="Downsampling with Seqtk",
        )

    class Output:
        """Output fields."""

        # Workflows do not have output fields.

    def run(self, inputs, outputs):
        """Run the workflow."""

        if not inputs.cdna_index and inputs.assay_type == "auto":
            self.error(
                "The input cDNA index file is necessary for 'Detect automatically' "
                "assay type."
            )

        input_bbduk = {
            "reads": inputs.reads,
            "min_length": inputs.preprocessing.min_length,
            "reference": {
                "sequences": inputs.preprocessing.adapters or [],
                "literal_sequences": inputs.preprocessing.custom_adapter_sequences,
            },
            "processing": {
                "kmer_length": inputs.preprocessing.kmer_length,
                "hamming_distance": inputs.preprocessing.hamming_distance,
            },
            "operations": {
                "quality_trim": "r",
                "trim_quality": inputs.preprocessing.trim_quality,
                "quality_encoding_offset": inputs.preprocessing.quality_encoding_offset,
                "ignore_bad_quality": inputs.preprocessing.ignore_bad_quality,
                "maxns": inputs.preprocessing.maxns,
            },
        }

        if (
            inputs.preprocessing.adapters
            or inputs.preprocessing.custom_adapter_sequences
        ):
            input_bbduk["operations"]["k_trim"] = "r"
        else:
            input_bbduk["operations"]["k_trim"] = "f"

        if (
            inputs.preprocessing.adapters
            or inputs.preprocessing.custom_adapter_sequences
        ):
            input_bbduk["operations"]["min_k"] = inputs.preprocessing.min_k
        else:
            input_bbduk["operations"]["min_k"] = -1

        if inputs.reads.type.startswith("data:reads:fastq:single:"):
            slug_bbduk = "bbduk-single"
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            input_bbduk["operations"]["trim_pairs_evenly"] = True
            input_bbduk["operations"]["trim_by_overlap"] = True

            slug_bbduk = "bbduk-paired"
        else:
            self.error("Wrong reads input type.")

        preprocessing = Data.create(
            process=BioProcess.get_latest(slug=slug_bbduk),
            input=input_bbduk,
            name=f"Trimmed ({inputs.reads.name})",
        )

        input_star = {
            "reads": preprocessing,
            "genome": inputs.genome,
            "unstranded": inputs.alignment.unstranded,
            "noncannonical": inputs.alignment.noncannonical,
            "gene_counts": True,
            "detect_chimeric": {
                "chimeric": inputs.alignment.chimeric_reads.chimeric,
                "chim_segment_min": inputs.alignment.chimeric_reads.chim_segment_min,
                "chim_out_type": inputs.alignment.chimeric_reads.chim_out_type,
                "chim_junction_overhang_min": inputs.alignment.chimeric_reads.chim_junction_overhang_min,
                "chim_score_drop_max": inputs.alignment.chimeric_reads.chim_score_drop_max,
                "chim_score_junction_non_gtag": inputs.alignment.chimeric_reads.chim_score_junction_non_gtag,
                "chim_score_separation": inputs.alignment.chimeric_reads.chim_score_separation,
                "chim_segment_read_gap_max": inputs.alignment.chimeric_reads.chim_segment_read_gap_max,
                "chim_multimap_nmax": inputs.alignment.chimeric_reads.chim_multimap_nmax,
                "pe_overlap_n_bases_min": inputs.alignment.chimeric_reads.pe_overlap_n_bases_min,
                "align_spliced_mate_map_lmin_over_lmate": inputs.alignment.chimeric_reads.align_spliced_mate_map_lmin_over_lmate,
                "align_sj_stitch_mismatch_n_max": inputs.alignment.chimeric_reads.align_sj_stitch_mismatch_n_max,
            },
            "t_coordinates": {
                "quant_mode": inputs.alignment.transcript_output.quant_mode,
                "single_end": inputs.alignment.transcript_output.single_end,
            },
            "filtering": {
                "out_filter_type": inputs.alignment.filtering_options.out_filter_type,
            },
            "alignment": {
                "align_end_alignment": inputs.alignment.alignment_options.align_end_alignment
            },
            "two_pass_mapping": {
                "two_pass_mode": inputs.alignment.two_pass_mapping.two_pass_mode
            },
            "output_options": {
                "out_unmapped": inputs.alignment.output_options.out_unmapped,
                "out_sam_attributes": inputs.alignment.output_options.out_sam_attributes,
            },
        }
        if inputs.alignment.filtering_options.out_multimap_max:
            input_star["filtering"][
                "out_multimap_max"
            ] = inputs.alignment.filtering_options.out_multimap_max

        if inputs.alignment.chimeric_reads.chimeric:
            input_star["filtering"]["out_multimap_max"] = 50

        if inputs.alignment.filtering_options.out_mismatch_max:
            input_star["filtering"][
                "out_mismatch_max"
            ] = inputs.alignment.filtering_options.out_mismatch_max

        if inputs.alignment.filtering_options.out_mismatch_nl_max:
            input_star["filtering"][
                "out_mismatch_nl_max"
            ] = inputs.alignment.filtering_options.out_mismatch_nl_max

        if inputs.alignment.filtering_options.out_score_min:
            input_star["filtering"][
                "out_score_min"
            ] = inputs.alignment.filtering_options.out_score_min

        if inputs.alignment.filtering_options.out_mismatch_nrl_max:
            input_star["filtering"][
                "out_mismatch_nrl_max"
            ] = inputs.alignment.filtering_options.out_mismatch_nrl_max

        if inputs.alignment.alignment_options.align_overhang_min:
            input_star["alignment"][
                "align_overhang_min"
            ] = inputs.alignment.alignment_options.align_overhang_min

        if inputs.alignment.alignment_options.align_sjdb_overhang_min:
            input_star["alignment"][
                "align_sjdb_overhang_min"
            ] = inputs.alignment.alignment_options.align_sjdb_overhang_min

        if inputs.alignment.alignment_options.align_intron_size_min:
            input_star["alignment"][
                "align_intron_size_min"
            ] = inputs.alignment.alignment_options.align_intron_size_min

        if inputs.alignment.alignment_options.align_intron_size_max:
            input_star["alignment"][
                "align_intron_size_max"
            ] = inputs.alignment.alignment_options.align_intron_size_max

        if inputs.alignment.alignment_options.align_gap_max:
            input_star["alignment"][
                "align_gap_max"
            ] = inputs.alignment.alignment_options.align_gap_max

        if inputs.alignment.output_options.out_rg_line:
            input_star["output_options"][
                "out_rg_line"
            ] = inputs.alignment.output_options.out_rg_line

        alignment = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input=input_star,
            name=f"Aligned ({inputs.reads.name})",
        )

        input_quant = {
            "aligned_reads": alignment,
            "n_reads": inputs.quantification.n_reads,
            "assay_type": inputs.assay_type,
            "annotation": inputs.annotation,
        }

        if inputs.cdna_index:
            input_quant["cdna_index"] = inputs.cdna_index

        counts = Data.create(
            process=BioProcess.get_latest(slug="star-quantification"),
            input=input_quant,
            name=f"Quantified ({inputs.reads.name})",
        )

        input_seqtk = {
            "reads": preprocessing,
            "n_reads": inputs.downsampling.n_reads,
            "advanced": {
                "seed": inputs.downsampling.advanced.seed,
                "fraction": inputs.downsampling.advanced.fraction,
                "two_pass": inputs.downsampling.advanced.two_pass,
            },
        }

        if inputs.reads.type.startswith("data:reads:fastq:single:"):
            slug_seqtk = "seqtk-sample-single"
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            slug_seqtk = "seqtk-sample-paired"
        else:
            self.error("Wrong reads input type.")

        downsampling = Data.create(
            process=BioProcess.get_latest(slug=slug_seqtk),
            input=input_seqtk,
            name=f"Subsampled ({inputs.reads.name})",
        )

        alignment_qc_rrna = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input={
                "reads": downsampling,
                "genome": inputs.rrna_reference,
            },
            name=f"rRNA aligned ({inputs.reads.name})",
        )
        alignment_qc_globin = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input={
                "reads": downsampling,
                "genome": inputs.globin_reference,
            },
            name=f"Globin aligned ({inputs.reads.name})",
        )
        idxstats = Data.create(
            process=BioProcess.get_latest(slug="samtools-idxstats"),
            input={
                "alignment": alignment,
            },
            name=f"Alignment summary ({inputs.reads.name})",
        )
        alignment_downsampled = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input={
                "reads": downsampling,
                "genome": inputs.genome,
            },
            name=f"Aligned subset ({inputs.reads.name})",
        )

        input_qorts = {
            "alignment": alignment_downsampled,
            "annotation": inputs.annotation,
            "options": {
                "stranded": inputs.assay_type,
            },
        }

        if inputs.cdna_index:
            input_qorts["options"]["cdna_index"] = inputs.cdna_index

        qorts = Data.create(
            process=BioProcess.get_latest(slug="qorts-qc"),
            input=input_qorts,
            name=f"QoRTs QC report ({inputs.reads.name})",
        )

        input_multiqc = {
            "data": [
                inputs.reads,
                preprocessing,
                alignment,
                counts,
                downsampling,
                alignment_qc_rrna,
                alignment_qc_globin,
                idxstats,
                qorts,
            ],
            "advanced": {
                "dirs": True,
                "config": True,
                "qc_evaluation": inputs.qc_evaluation,
            },
        }

        Data.create(process=BioProcess.get_latest(slug="multiqc"), input=input_multiqc)

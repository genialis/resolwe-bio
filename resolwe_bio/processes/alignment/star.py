"""Align reads with STAR aligner."""

import gzip
import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)

SPECIES = [
    "Caenorhabditis elegans",
    "Cricetulus griseus",
    "Dictyostelium discoideum",
    "Dictyostelium purpureum",
    "Drosophila melanogaster",
    "Homo sapiens",
    "Macaca mulatta",
    "Mus musculus",
    "Rattus norvegicus",
    "other",
]


def get_fastq_name(fastq_path):
    """Get the name of the FASTQ file."""
    fastq_file = fastq_path.name
    assert fastq_file.endswith(".fastq.gz")
    return fastq_file[:-9]


class AlignmentStar(Process):
    """Align reads with STAR aligner.

    Spliced Transcripts Alignment to a Reference (STAR) software is
    based on an alignment algorithm that uses sequential maximum
    mappable seed search in uncompressed suffix arrays followed by seed
    clustering and stitching procedure. In addition to unbiased de novo
    detection of canonical junctions, STAR can discover non-canonical
    splices and chimeric (fusion) transcripts, and is also capable of
    mapping full-length RNA sequences. More information can be found in
    the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
    and in the [original paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/).
    The current version of STAR is 2.7.10b.
    """

    slug = "alignment-star"
    name = "STAR"
    process_type = "data:alignment:bam:star"
    version = "5.1.1"
    category = "Align"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.2.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
        },
    }
    data_name = "{{ reads|name|default('?') }}"

    class Input:
        """Input fields to process AlignmentStar."""

        reads = DataField("reads:fastq", label="Input reads (FASTQ)")

        genome = DataField(
            "index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
        )

        annotation = DataField(
            "annotation",
            label="Annotation file (GTF/GFF3)",
            required=False,
            description="Insert known annotations into genome indices at the mapping stage.",
        )

        unstranded = BooleanField(
            label="The data is unstranded [--outSAMstrandField intronMotif]",
            default=False,
            description="For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced "
            "alignments with XS strand attribute, which STAR will generate with "
            "--outSAMstrandField intronMotif option. As required, the XS strand attribute will be "
            "generated for all alignments that contain splice junctions. The spliced alignments "
            "that have undefined strand (i.e. containing only non-canonical unannotated "
            "junctions) will be suppressed. If you have stranded RNA-seq data, you do not need to "
            "use any specific STAR options. Instead, you need to run Cufflinks with the library "
            "option --library-type options. For example, cufflinks --library-type fr-firststrand "
            "should be used for the standard dUTP protocol, including Illumina's stranded "
            "Tru-Seq. This option has to be used only for Cufflinks runs and not for STAR runs.",
        )

        noncannonical = BooleanField(
            label="Remove non-canonical junctions (Cufflinks compatibility)",
            default=False,
            description="It is recommended to remove the non-canonical junctions for Cufflinks "
            "runs using --outFilterIntronMotifs RemoveNoncanonical.",
        )

        gene_counts = BooleanField(
            label="Gene count [--quantMode GeneCounts]",
            description="With this option set to True STAR will count the number of reads per gene "
            "while mapping. A read is counted if it overlaps (1nt or more) one and only one "
            "gene. Both ends of the paired-end read are checked for overlaps. The counts coincide "
            "with those produced by htseq-count with default parameters.",
            default=False,
        )

        class AnnotationOptions:
            """Annotation file options."""

            feature_exon = StringField(
                label="Feature type [--sjdbGTFfeatureExon]",
                default="exon",
                description="Feature type in GTF file to be used as exons for building "
                "transcripts.",
            )

            sjdb_overhang = IntegerField(
                label="Junction length [--sjdbOverhang]",
                default=100,
                description="This parameter specifies the length of the genomic sequence around "
                "the annotated junction to be used in constructing the splice junction database. "
                "Ideally, this length should be equal to the ReadLength-1, where ReadLength is "
                "the length of the reads. For instance, for Illumina 2x100b paired-end reads, the "
                "ideal value is 100-1=99. In the case of reads of varying length, the ideal value "
                "is max(ReadLength)-1. In most cases, the default value of 100 will work as well "
                "as the ideal value.",
            )

        class ChimericReadsOptions:
            """Chimeric reads options."""

            chimeric = BooleanField(
                label="Detect chimeric and circular alignments [--chimOutType SeparateSAMold]",
                default=False,
                description="To switch on detection of chimeric (fusion) alignments (in addition "
                "to normal mapping), --chimSegmentMin should be set to a positive value. Each "
                "chimeric alignment consists of two segments.Each segment is non-chimeric on "
                "its own, but the segments are chimeric to each other (i.e. the segments belong "
                "to different chromosomes, or different strands, or are far from each other). "
                "Both segments may contain splice junctions, and one of the segments may contain "
                "portions of both mates. --chimSegmentMin parameter controls the minimum mapped "
                "length of the two segments that is allowed. For example, if you have 2x75 reads "
                "and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome "
                "and 20b on the other will be output, while 135 + 15 won't be.",
            )

            chim_segment_min = IntegerField(
                label="Minimum length of chimeric segment [--chimSegmentMin]",
                default=20,
                disabled="!detect_chimeric.chimeric",
            )

        class TranscriptOutputOptions:
            """Transcript coordinate output options."""

            quant_mode = BooleanField(
                label="Output in transcript coordinates [--quantMode TranscriptomeSAM]",
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
                label="Minumum alignment score [--outFilterScoreMin]",
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
            """Alignment options."""

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
                required=False,
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
            )

        class TwoPassOptions:
            """Two-pass mapping options."""

            two_pass_mode = BooleanField(
                label="Use two pass mode [--twopassMode]",
                default=False,
                description="Use two-pass maping instead of first-pass only. In two-pass mode we "
                "first perform first-pass mapping, extract junctions, insert them into genome "
                "index, and re-map all reads in the second mapping pass.",
            )

        class OutputOptions:
            """Output options."""

            out_unmapped = BooleanField(
                label="Output unmapped reads (SAM) [--outSAMunmapped Within]",
                default=False,
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
                "quoted. Comma separated RG lines correspons to different (comma separated) input "
                "files in –readFilesIn. Commas have to be surrounded by spaces, e.g. "
                "–outSAMattrRGline ID:xxx , ID:zzz ”DS:z z” , ID:yyy DS:yyyy.",
            )

        class Limits:
            """Limits."""

            limit_buffer_size = ListField(
                IntegerField(),
                label="Buffer size [--limitIObufferSize]",
                default=[30000000, 50000000],
                description="Maximum available buffers size (bytes) for input/output, per thread. "
                "Parameter requires two numbers - separate sizes for input and output buffers.",
            )

            limit_sam_records = IntegerField(
                label="Maximum size of the SAM record [--limitOutSAMoneReadBytes]",
                default=100000,
                description="Maximum size of the SAM record (bytes) for one read. Recommended "
                "value: >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax.",
            )

            limit_junction_reads = IntegerField(
                label="Maximum number of junctions [--limitOutSJoneRead]",
                default=1000,
                description="Maximum number of junctions for one read (including all "
                "multi-mappers).",
            )

            limit_collapsed_junctions = IntegerField(
                label="Maximum number of collapsed junctions [--limitOutSJcollapsed]",
                default=1000000,
            )

            limit_inserted_junctions = IntegerField(
                label="Maximum number of junction to be inserted [--limitSjdbInsertNsj]",
                default=1000000,
                description="Maximum number of junction to be inserted to the genome on the fly  "
                "at the mapping stage, including those from annotations and those detected in the "
                "1st step of the 2-pass run.",
            )

        annotation_options = GroupField(
            AnnotationOptions, label="Annotation file options", hidden="!annotation"
        )
        detect_chimeric = GroupField(
            ChimericReadsOptions, label="Chimeric and circular alignments"
        )
        t_coordinates = GroupField(
            TranscriptOutputOptions, label="Transcript coordinates output"
        )
        filtering = GroupField(FilteringOptions, label="Output Filtering")
        alignment = GroupField(AlignmentOptions, label="Alignment and Seeding")
        two_pass_mapping = GroupField(TwoPassOptions, label="Two-pass mapping")
        output_options = GroupField(OutputOptions, label="Output options")
        limits = GroupField(Limits, label="Limits")

    class Output:
        """Output fields to process AlignmentStar."""

        bam = FileField(label="Alignment file")
        bai = FileField(label="BAM file index")
        unmapped_1 = FileField(label="Unmapped reads (mate 1)", required=False)
        unmapped_2 = FileField(label="Unmapped reads (mate 2)", required=False)
        sj = FileField(label="Splice junctions")
        chimeric = FileField(label="Chimeric alignments", required=False)
        alignment_transcriptome = FileField(
            label="Alignment (transcriptome coordinates)", required=False
        )
        gene_counts = FileField(label="Gene counts", required=False)
        stats = FileField(label="Statistics")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""

        reads_species = inputs.reads.entity.annotations["general.species"]
        if reads_species is not None and reads_species != inputs.genome.output.species:
            self.warning(
                f"Species of reads ({inputs.reads.entity.annotations['general.species']}) "
                f"and genome ({inputs.genome.output.species}) do not match."
            )

        elif inputs.genome.output.species in SPECIES and reads_species is None:
            inputs.reads.entity.annotations["general.species"] = (
                inputs.genome.output.species
            )
            self.info("Sample species was automatically annotated to match the genome.")

        mate1_name = get_fastq_name(Path(inputs.reads.output.fastq[0].path))
        mate_1 = [fastq.path for fastq in inputs.reads.output.fastq]

        if inputs.reads.type.startswith("data:reads:fastq:paired:"):
            mate2_name = get_fastq_name(Path(inputs.reads.output.fastq2[0].path))
            mate_2 = [fastq.path for fastq in inputs.reads.output.fastq2]

        self.progress(0.05)

        star_params = [
            "--runThreadN",
            self.requirements.resources.cores,
            "--genomeDir",
            inputs.genome.output.index.path,
            "--outReadsUnmapped",
            "Fastx",
            "--limitIObufferSize",
            inputs.limits.limit_buffer_size[0],
            inputs.limits.limit_buffer_size[1],
            "--limitOutSAMoneReadBytes",
            inputs.limits.limit_sam_records,
            "--limitOutSJoneRead",
            inputs.limits.limit_junction_reads,
            "--limitOutSJcollapsed",
            inputs.limits.limit_collapsed_junctions,
            "--limitSjdbInsertNsj",
            inputs.limits.limit_inserted_junctions,
            "--outFilterType",
            inputs.filtering.out_filter_type,
            "--outSAMtype",
            "BAM",
            "Unsorted",
        ]

        if inputs.reads.type.startswith("data:reads:fastq:single:"):
            star_params.extend(
                ["--readFilesIn", ",".join(mate_1), "--readFilesCommand", "zcat"]
            )
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            star_params.extend(
                [
                    "--readFilesIn",
                    ",".join(mate_1),
                    ",".join(mate_2),
                    "--readFilesCommand",
                    "zcat",
                ]
            )
        else:
            self.error("Wrong reads input type.")

        if inputs.annotation:
            star_params.extend(
                [
                    "--sjdbGTFfile",
                    inputs.annotation.output.annot.path,
                    "--sjdbOverhang",
                    inputs.annotation_options.sjdb_overhang,
                    "--sjdbGTFfeatureExon",
                    inputs.annotation_options.feature_exon,
                ]
            )

            if inputs.annotation.type.startswith("data:annotation:gff3:"):
                star_params.extend(["--sjdbGTFtagExonParentTranscript", "Parent"])

        if inputs.unstranded:
            star_params.extend(["--outSAMstrandField", "intronMotif"])

        if inputs.noncannonical:
            star_params.extend(["--outFilterIntronMotifs", "RemoveNoncanonical"])

        if inputs.detect_chimeric.chimeric:
            star_params.extend(
                [
                    "--chimOutType",
                    "SeparateSAMold",
                    "--chimSegmentMin",
                    inputs.detect_chimeric.chim_segment_min,
                ]
            )

        gene_segments = Path(inputs.genome.output.index.path) / "geneInfo.tab"
        if inputs.t_coordinates.quant_mode:
            if not gene_segments.is_file() and not inputs.annotation:
                self.error(
                    "Output in transcript coordinates requires genome annotation file."
                )

            if inputs.gene_counts:
                star_params.extend(["--quantMode", "TranscriptomeSAM", "GeneCounts"])
            else:
                star_params.extend(["--quantMode", "TranscriptomeSAM"])

            if inputs.t_coordinates.single_end:
                star_params.extend(["--quantTranscriptomeBan", "Singleend"])

        elif inputs.gene_counts:
            if not gene_segments.is_file() and not inputs.annotation:
                self.error(
                    "Counting the number of reads per gene requires a genome "
                    "annotation file."
                )

            star_params.extend(["--quantMode", "GeneCounts"])

        if inputs.filtering.out_multimap_max:
            star_params.extend(
                ["--outFilterMultimapNmax", inputs.filtering.out_multimap_max]
            )

        if inputs.filtering.out_mismatch_max:
            star_params.extend(
                ["--outFilterMismatchNmax", inputs.filtering.out_mismatch_max]
            )

        if inputs.filtering.out_mismatch_nl_max:
            star_params.extend(
                ["--outFilterMismatchNoverLmax", inputs.filtering.out_mismatch_nl_max]
            )

        if inputs.filtering.out_score_min:
            star_params.extend(["--outFilterScoreMin", inputs.filtering.out_score_min])

        if inputs.filtering.out_mismatch_nrl_max:
            star_params.extend(
                [
                    "--outFilterMismatchNoverReadLmax",
                    inputs.filtering.out_mismatch_nrl_max,
                ]
            )

        if inputs.alignment.align_overhang_min:
            star_params.extend(
                ["--alignSJoverhangMin", inputs.alignment.align_overhang_min]
            )

        if inputs.alignment.align_sjdb_overhang_min:
            star_params.extend(
                ["--alignSJDBoverhangMin", inputs.alignment.align_sjdb_overhang_min]
            )

        if inputs.alignment.align_intron_size_min:
            star_params.extend(
                ["--alignIntronMin", inputs.alignment.align_intron_size_min]
            )

        if inputs.alignment.align_intron_size_max:
            star_params.extend(
                ["--alignIntronMax", inputs.alignment.align_intron_size_max]
            )

        if inputs.alignment.align_gap_max:
            star_params.extend(["--alignMatesGapMax", inputs.alignment.align_gap_max])

        if inputs.alignment.align_end_alignment:
            star_params.extend(
                ["--alignMatesGapMax", inputs.alignment.align_end_alignment]
            )

        if inputs.two_pass_mapping.two_pass_mode:
            star_params.extend(["--twopassMode", "Basic"])

        if inputs.output_options.out_unmapped:
            star_params.extend(["--outSAMunmapped", "Within"])

        if inputs.output_options.out_sam_attributes:
            # Create a list from string of out_sam_attributes to avoid unknown/unimplemented
            # SAM attrribute error due to Plumbum command passing problems.
            attributes = inputs.output_options.out_sam_attributes.split(" ")
            star_params.extend(["--outSAMattributes", attributes])

        if inputs.output_options.out_rg_line:
            star_params.extend(
                ["--outSAMattrRGline", inputs.output_options.out_rg_line]
            )
        elif len(mate_1) > 1:
            read_groups = [
                f"ID:{Path(file_path).name[:-9].replace(' ', '_')} SM:sample1"
                for file_path in mate_1
            ]
            star_params.append("--outSAMattrRGline " + " , ".join(read_groups))

        self.progress(0.1)

        return_code, _, _ = Cmd["STAR"][star_params] & TEE(retcode=None)
        log_file = Path("Log.out")
        # Log contains useful information for debugging.
        if log_file.is_file():
            with open(log_file, "r") as log:
                print(log.read())
        if return_code:
            self.error("Reads alignment failed.")

        self.progress(0.7)

        star_unmapped_r1 = Path("Unmapped.out.mate1")
        if star_unmapped_r1.is_file():
            unmapped_out_1 = f"{mate1_name}_unmapped.out.mate1.fastq"
            star_unmapped_r1.rename(unmapped_out_1)

            return_code, _, _ = Cmd["pigz"][unmapped_out_1] & TEE(retcode=None)
            if return_code:
                self.error("Compression of unmapped mate 1 reads failed.")

            outputs.unmapped_1 = f"{unmapped_out_1}.gz"

        star_unmapped_r2 = Path("Unmapped.out.mate2")
        if (
            inputs.reads.type.startswith("data:reads:fastq:paired:")
            and star_unmapped_r2.is_file()
        ):
            unmapped_out_2 = f"{mate2_name}_unmapped.out.mate2.fastq"
            star_unmapped_r2.rename(unmapped_out_2)

            return_code, _, _ = Cmd["pigz"][unmapped_out_2] & TEE(retcode=None)
            if return_code:
                self.error("Compression of unmapped mate 2 reads failed.")
            outputs.unmapped_2 = f"{unmapped_out_2}.gz"

        self.progress(0.8)

        out_bam = f"{mate1_name}.bam"
        out_bai = f"{out_bam}.bai"

        sort_params = [
            "Aligned.out.bam",
            "-o",
            out_bam,
            "-@",
            self.requirements.resources.cores,
        ]
        return_code, _, _ = Cmd["samtools"]["sort"][sort_params] & TEE(retcode=None)
        if return_code:
            self.error("Samtools sort command failed.")

        outputs.bam = out_bam

        return_code, _, _ = Cmd["samtools"]["index"][out_bam, out_bai] & TEE(
            retcode=None
        )
        if return_code:
            self.error("Samtools index command failed.")

        outputs.bai = out_bai

        self.progress(0.9)

        if inputs.detect_chimeric.chimeric:
            out_chimeric = f"{mate1_name}_chimeric.out.sam"
            Path("Chimeric.out.sam").rename(out_chimeric)
            outputs.chimeric = out_chimeric

        if inputs.t_coordinates.quant_mode:
            out_transcriptome = f"{mate1_name}_aligned.toTranscriptome.out.bam"
            Path("Aligned.toTranscriptome.out.bam").rename(out_transcriptome)
            outputs.alignment_transcriptome = out_transcriptome

        if inputs.gene_counts:
            out_counts = f"{mate1_name}_ReadsPerGene.out.tab.gz"
            with open(file="ReadsPerGene.out.tab", mode="rb") as f_in:
                with gzip.open(filename=out_counts, mode="wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            outputs.gene_counts = out_counts

        out_stats = f"{mate1_name}_stats.txt"
        Path("Log.final.out").rename(out_stats)
        outputs.stats = out_stats

        out_sj = f"{mate1_name}_SJ.out.tab"
        Path("SJ.out.tab").rename(out_sj)
        outputs.sj = out_sj

        outputs.species = inputs.genome.output.species
        outputs.build = inputs.genome.output.build

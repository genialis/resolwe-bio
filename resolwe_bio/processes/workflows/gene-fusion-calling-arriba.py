"""Gene Fusion Calling with Arriba."""

from resolwe.process import (
    BooleanField,
    Data,
    DataField,
    FloatField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)
from resolwe.process.models import Process as BioProcess


class WorkflowGeneFusionCallingArriba(Process):
    """Gene Fusion Calling with Arriba.

    This workflow is designed to detect gene fusions from RNA-Seq data using the Arriba tool.
    The workflow accommodates two types of input:
    - raw RNA-Seq reads in FASTQ format or
    - pre-aligned BAM files generated by the STAR aligner.

    1. **Input Requirements:**
        - Users can provide either FASTQ files or BAM files as input, but not both simultaneously.
        - If FASTQ files are provided, the workflow will align these reads using the STAR aligner.
        - The workflow does not perform quality trimming (e.g., BBDuk) on FASTQ inputs, assuming
          that the reads are already preprocessed.
        - If BAM files are provided, the workflow checks if they were generated by STAR aligner
          and confirms that the appropriate alignment settings were used (handled by Arriba).
    """

    slug = "gene-fusion-calling-arriba"
    name = "Gene Fusion Calling with Arriba"
    version = "1.1.1"
    process_type = "data:workflow:genefusions:arriba"
    category = "Pipeline"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default(bam|name|default('?')) }}"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Reads (FASTQ)",
            description="Preprocessed paired end reads in FASTQ format.",
            required=False,
            disabled="bam",
        )

        bam = DataField(
            data_type="alignment:bam",
            label="Input BAM file from STAR ran with parameters suggested by Arriba.",
            required=False,
            disabled="reads",
        )

        gtf = DataField(
            data_type="annotation:gtf",
            label="GTF file",
            description="Annotation file in GTF format.",
        )

        genome = DataField(
            data_type="seq:nucleotide",
            label="Genome file",
            description="Genome file in FASTA format.",
        )

        star_index = DataField(
            data_type="index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
            required=False,
            hidden="!reads",
        )

        blacklist_file = DataField(
            data_type="file",
            label="Blacklist file",
            description="Arriba blacklist file.",
            required=False,
        )

        known_fusions_file = DataField(
            data_type="file",
            label="Known fusions file",
            description="Arriba known fusions file.",
            required=False,
        )

        class ChimericReadsOptions:
            """Chimeric reads options."""

            chimeric = BooleanField(
                label="Detect chimeric and circular alignments [--chimOutType]",
                default=True,
                description="Detect chimeric and circular alignments. If set to True, "
                "STAR parameters are set to values recommended for maximum sensitivity "
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
                description="Type of chimeric output produced by STAR.",
            )

            chim_segment_min = IntegerField(
                label="Minimum length of chimeric segment [--chimSegmentMin]",
                default=10,
                description="Minimum length of chimeric segment to be reported.",
            )

            chim_junction_overhang_min = IntegerField(
                label="Minimum overhang for chimeric junctions [--chimJunctionOverhangMin]",
                default=10,
                description="Minimum overhang for chimeric junctions.",
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
            )

            chim_score_separation = IntegerField(
                label="Minimum score separation between chimeric alignments [--chimScoreSeparation]",
                default=1,
                description="Minimum score separation between the best chimeric alignment "
                "and the next best.",
            )

            chim_segment_read_gap_max = IntegerField(
                label="Maximum read gap for chimeric segments [--chimSegmentReadGapMax]",
                default=3,
                description="Maximum gap in the read sequence between chimeric segments.",
            )

            chim_multimap_nmax = IntegerField(
                label="Maximum number of chimeric multi-mapping segments [--chimMultimapNmax]",
                default=50,
                description="Maximum number of chimeric alignments reported for multimapped reads.",
            )

            pe_overlap_n_bases_min = IntegerField(
                label="Minimum number of overlapping bases in PE [--peOverlapNbasesMin]",
                default=10,
                description="Minimum number of overlapping bases to trigger merging mate pair.",
            )

            align_spliced_mate_map_lmin_over_lmate = FloatField(
                label="Minimum mapped length of spliced mates over mate "
                "length [--alignSplicedMateMapLminOverLmate]",
                default=0.5,
                description="Minimum mapped length of the spliced mate, divided by the mate length.",
            )

            align_sj_stitch_mismatch_n_max = ListField(
                IntegerField(),
                label="Maximum number of mismatches for stitching SJ [--alignSJstitchMismatchNmax]",
                default=[5, -1, 5, 5],
                description="Maximum number of mismatches that are allowed in the seed region for "
                "stitching splice junctions. Enter four integers.",
            )

            out_multimap_max = IntegerField(
                label="Maximum number of loci [--outFilterMultimapNmax]",
                required=False,
                description="Maximum number of loci the read is allowed to map to. Alignments "
                "(all of them) will be output only if the read maps to no more loci than this "
                "value. Otherwise no alignments will be output, and the read will be counted as "
                "'mapped to too many loci' (default: 10).",
            )

        star_chimeric = GroupField(
            ChimericReadsOptions,
            label="Alignment with STAR",
        )

    class Output:
        """Output fields."""

    def run(self, inputs, outputs):
        """Run the workflow."""

        if inputs.reads and inputs.bam:
            self.error("Select either FASTQ or BAM file, not both.")

        if not inputs.reads and not inputs.bam:
            self.error("You must provide either FASTQ or BAM input.")

        if inputs.bam:
            bam = inputs.bam

        elif inputs.reads:
            input_star = {
                "reads": inputs.reads,
                "genome": inputs.star_index,
                "gene_counts": True,
                "detect_chimeric": {
                    "chimeric": inputs.star_chimeric.chimeric,
                    "chim_segment_min": inputs.star_chimeric.chim_segment_min,
                    "chim_out_type": inputs.star_chimeric.chim_out_type,
                    "chim_junction_overhang_min": inputs.star_chimeric.chim_junction_overhang_min,
                    "chim_score_drop_max": inputs.star_chimeric.chim_score_drop_max,
                    "chim_score_junction_non_gtag": inputs.star_chimeric.chim_score_junction_non_gtag,
                    "chim_score_separation": inputs.star_chimeric.chim_score_separation,
                    "chim_segment_read_gap_max": inputs.star_chimeric.chim_segment_read_gap_max,
                    "chim_multimap_nmax": inputs.star_chimeric.chim_multimap_nmax,
                    "pe_overlap_n_bases_min": inputs.star_chimeric.pe_overlap_n_bases_min,
                    "align_spliced_mate_map_lmin_over_lmate": inputs.star_chimeric.align_spliced_mate_map_lmin_over_lmate,
                    "align_sj_stitch_mismatch_n_max": inputs.star_chimeric.align_sj_stitch_mismatch_n_max,
                },
                "filtering": {
                    "out_multimap_max": inputs.star_chimeric.out_multimap_max,
                },
            }

            bam = Data.create(
                process=BioProcess.get_latest(slug="alignment-star"),
                input=input_star,
                name=f"Aligned ({inputs.reads.name})",
            )

        arriba_input = {
            "bam": bam.id,
            "gtf": inputs.gtf.id,
            "genome": inputs.genome.id,
        }

        if inputs.blacklist_file:
            arriba_input["blacklist_file"] = inputs.blacklist_file.id

        if inputs.known_fusions_file:
            arriba_input["known_fusions_file"] = inputs.known_fusions_file.id

        fusion_name = (
            f"Fusions ({inputs.bam.name})"
            if inputs.bam
            else f"Fusions ({inputs.reads.name})"
        )

        Data.create(
            process=BioProcess.get_latest(slug="arriba"),
            input=arriba_input,
            name=fusion_name,
        )

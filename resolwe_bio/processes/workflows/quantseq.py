"""3' mRNA-Seq pipeline comprised of QC, preprocessing, alignment and quantification."""

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


class WorkflowQuantSeq(Process):
    """3' mRNA-Seq pipeline.

    Reads are preprocessed by __BBDuk__ or __Cutadapt__ which removes adapters,
    trims reads for quality from the 3'-end, and discards reads that are too
    short after trimming. Preprocessed reads are aligned by __STAR__
    aligner. For read-count quantification, the __FeatureCounts__ tool
    is used. QoRTs QC and Samtools idxstats tools are used to report
    alignment QC metrics.

    QC steps include downsampling, QoRTs QC analysis and alignment of
    input reads to the rRNA/globin reference sequences. The reported
    alignment rate is used to assess the rRNA/globin sequence depletion
    rate.
    """

    slug = "workflow-quantseq"
    name = "QuantSeq workflow"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    version = "5.1.0"
    entity = {
        "type": "sample",
    }
    process_type = "data:workflow:quant:featurecounts"
    category = "Pipeline"

    class Input:
        """Input fields."""

        trimming_tool = StringField(
            label="Trimming tool",
            choices=[
                ("bbduk", "BBDuk"),
                ("cutadapt", "Cutadapt"),
            ],
            description="Select the trimming tool. If you select BBDuk then "
            "please provide adapter sequences in fasta file(s). If you select Cutadapt "
            "as a trimming tool, pre-determined adapter sequences will be removed.",
        )
        reads = DataField(
            data_type="reads:fastq",
            label="Input reads (FASTQ)",
            description="Reads in FASTQ file, single or paired end.",
        )
        genome = DataField(
            data_type="index:star",
            label="Indexed reference genome",
            description="Genome index prepared by STAR aligner indexing tool.",
        )
        adapters = ListField(
            inner=DataField(data_type="seq:nucleotide"),
            label="Adapters",
            required=False,
            hidden="trimming_tool != 'bbduk'",
            description="Provide a list of sequencing adapters files (.fasta) "
            "to be removed by BBDuk.",
        )
        annotation = DataField(
            data_type="annotation",
            label="Annotation",
            description="GTF and GFF3 annotation formats are supported.",
        )
        assay_type = StringField(
            label="Assay type",
            choices=[
                ("forward", "Strand-specific forward"),
                ("reverse", "Strand-specific reverse"),
            ],
            description="In strand-specific forward assay and single "
            "reads, the read has to be mapped to the same strand as the feature. "
            "For paired-end reads, the first read has to be on the same strand and "
            "the second read on the opposite strand. In strand-specific reverse "
            "assay these rules are reversed.",
            required=False,
        )
        rrna_reference = DataField(
            data_type="index:star",
            label="Indexed rRNA reference sequence",
            description="Reference sequence index prepared by STAR aligner indexing tool.",
            required=False,
        )
        globin_reference = DataField(
            data_type="index:star",
            label="Indexed Globin reference sequence",
            description="Reference sequence index prepared by STAR aligner indexing tool.",
            required=False,
        )

        class Preprocessing:
            """Preprocessing with BBDuk."""

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

        class Cutadapt:
            """Cutadapt filtering."""

            quality_cutoff = IntegerField(
                label="Reads quality cutoff",
                required=False,
                description="Trim low-quality bases from 3' end of each read before "
                "adapter removal. The use of this option will override the use "
                "of NextSeq/NovaSeq-specific trim option.",
            )

        class Downsampling:
            """Downsampling (Seqtk)."""

            n_reads = IntegerField(
                label="Number of reads",
                default=1000000,
                description="Number of reads to include in subsampling.",
            )

            class Advanced:
                """Advanced options for downsampling."""

                seed = IntegerField(
                    label="Number of reads",
                    default=11,
                    description="Using the same random seed makes reads subsampling reproducible "
                    "in different environments.",
                )
                fraction = FloatField(
                    label="Fraction",
                    required=False,
                    range=[0, 1.0],
                    description="Use the fraction of reads [0 - 1.0] from the orignal input file instead "
                    "of the absolute number of reads. If set, this will override the"
                    "'Number of reads' input parameter.",
                )
                two_pass = BooleanField(
                    label="2-pass mode",
                    default=False,
                    description="Enable two-pass mode when down-sampling. Two-pass mode is twice "
                    "as slow but with much reduced memory.",
                )

            advanced = GroupField(Advanced, label="Advanced options for downsampling")

        cutadapt = GroupField(
            Cutadapt, label="Cutadapt filtering", hidden="trimming_tool != 'cutadapt'"
        )
        downsampling = GroupField(
            Downsampling,
            label="Downsampling with Seqtk",
        )
        preprocessing = GroupField(
            Preprocessing,
            label="Preprocessing with BBDuk",
            hidden="trimming_tool != 'bbduk'",
        )

    class Output:
        """Output fields."""

    def run(self, inputs, outputs):
        """Run the workflow."""

        if inputs.trimming_tool == "bbduk" and not inputs.adapters:
            self.error(
                "Please provide fasta file of adapters, if you want to use BBDuk as a trimming tool."
            )

        if inputs.trimming_tool == "bbduk":
            input_preprocessing = {
                "reads": inputs.reads,
                "min_length": 20,
                "reference": {"sequences": inputs.adapters},
                "processing": {"kmer_length": 13},
                "operations": {
                    "k_trim": "r",
                    "min_k": 6,
                    "quality_trim": "r",
                    "trim_quality": 10,
                    "quality_encoding_offset": inputs.preprocessing.quality_encoding_offset,
                    "ignore_bad_quality": inputs.preprocessing.ignore_bad_quality,
                },
                "fastqc": {
                    "nogroup": True,
                },
            }

            if inputs.reads.type.startswith("data:reads:fastq:single:"):
                process_slug = "bbduk-single"
            elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
                process_slug = "bbduk-paired"
            else:
                self.error("Wrong reads input type was provided.")

        else:
            if inputs.reads.type.startswith("data:reads:fastq:single:"):
                input_preprocessing = {
                    "reads": inputs.reads,
                    "options": {
                        "quality_cutoff": inputs.cutadapt.quality_cutoff,
                    },
                }
                process_slug = "cutadapt-3prime-single"
            else:
                self.error(
                    "Only single-end reads are supported when Cutadapt is "
                    "selected as a trimming tool."
                )

        preprocessing = Data.create(
            process=BioProcess.get_latest(process_slug),
            input=input_preprocessing,
            name=f"Trimmed ({inputs.reads.name})",
        )

        input_star = {
            "reads": preprocessing,
            "genome": inputs.genome,
            "filtering": {
                "out_filter_type": "BySJout",
                "out_multimap_max": 20,
                "out_mismatch_max": 999,
                "out_mismatch_nl_max": 0.6,
            },
            "alignment": {
                "align_overhang_min": 8,
                "align_sjdb_overhang_min": 1,
                "align_intron_size_min": 20,
                "align_intron_size_max": 1000000,
                "align_gap_max": 1000000,
            },
        }

        if inputs.trimming_tool == "cutadapt":
            input = {
                "unstranded": True,
                "output_options": {"out_sam_attributes": "NH HI NM MD"},
            }
            input_star.update(input)

        alignment = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input=input_star,
            name=f"Aligned ({inputs.reads.name})",
        )

        input_featurecounts = {
            "aligned_reads": alignment,
            "normalization_type": "CPM",
            "assay_type": inputs.assay_type,
            "annotation": inputs.annotation,
        }

        if inputs.trimming_tool == "cutadapt":
            input = {"assay_type": "forward"}
            input_featurecounts.update(input)

        quantification = Data.create(
            process=BioProcess.get_latest(slug="feature_counts"),
            input=input_featurecounts,
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
            process_slug = "seqtk-sample-single"
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            process_slug = "seqtk-sample-paired"
        else:
            self.error("Wrong reads input type was provided.")

        downsampling = Data.create(
            process=BioProcess.get_latest(slug=process_slug),
            input=input_seqtk,
            name=f"Subsampled ({inputs.reads.name})",
        )

        idxstats = Data.create(
            process=BioProcess.get_latest(slug="samtools-idxstats"),
            input={
                "alignment": alignment,
            },
            name=f"Alignment summary ({inputs.reads.name})",
        )

        alignment_qorts = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input={
                "reads": downsampling,
                "genome": inputs.genome,
            },
            name=f"Aligned subset ({inputs.reads.name})",
        )
        input_qorts = {
            "alignment": alignment_qorts,
            "annotation": inputs.annotation,
            "options": {
                "stranded": inputs.assay_type,
            },
        }

        if inputs.trimming_tool == "cutadapt":
            input = {"options": {"stranded": "forward"}}
            input_qorts.update(input)

        qorts = Data.create(
            process=BioProcess.get_latest(slug="qorts-qc"),
            input=input_qorts,
            name=f"QoRTs QC report ({inputs.reads.name})",
        )

        multiqc_inputs = [
            inputs.reads,
            preprocessing,
            alignment,
            downsampling,
            quantification,
            idxstats,
            qorts,
        ]

        if inputs.rrna_reference and inputs.globin_reference:
            alignment_qc_rrna = Data.create(
                process=BioProcess.get_latest(slug="alignment-star"),
                input={
                    "reads": downsampling,
                    "genome": inputs.rrna_reference,
                },
            )
            alignment_qc_globin = Data.create(
                process=BioProcess.get_latest(slug="alignment-star"),
                input={
                    "reads": downsampling,
                    "genome": inputs.globin_reference,
                },
            )

            multiqc_inputs.extend([alignment_qc_rrna, alignment_qc_globin])

        inputs_multiqc = {"data": multiqc_inputs}

        Data.create(process=BioProcess.get_latest(slug="multiqc"), input=inputs_multiqc)

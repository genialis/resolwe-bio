"""Alignment-free RNA-Seq pipeline."""

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


class WorkflowBbdukSalmonQc(Process):
    """Alignment-free RNA-Seq pipeline.

    Salmon tool and tximport package are used in quantification step to
    produce gene-level abundance estimates.

    rRNA and globin-sequence contamination rate in the sample is
    determined using STAR aligner. Quality-trimmed reads are down-sampled
    (using Seqtk tool) and aligned to the genome, rRNA and globin
    reference sequences. The rRNA and globin-sequence alignment rates
    indicate the percentage of the reads in the sample that are of
    rRNA and globin origin, respectively. Alignment of down-sampled data
    to a whole genome reference sequence is used to produce an alignment
    file suitable for Samtools and QoRTs QC analysis.

    Per-sample analysis results and QC data is summarized by the MultiQC
    tool.
    """

    slug = "workflow-bbduk-salmon-qc"
    name = "BBDuk - Salmon - QC"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    entity = {
        "type": "sample",
    }
    version = "4.4.0"
    process_type = "data:workflow:rnaseq:salmon"
    category = "Pipeline"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Select sample(s) (FASTQ)",
            description="Reads in FASTQ file, single or paired end.",
        )
        salmon_index = DataField(
            data_type="index:salmon",
            label="Salmon index",
            description="Transcriptome index file created using the Salmon indexing tool.",
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

        class Preprocessing:
            """Preprocessing with BBDuk."""

            adapters = ListField(
                inner=DataField(data_type="seq:nucleotide"),
                label="Adapters",
                description="FASTA file(s) with adapters.",
                required=False,
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
                label="K-mer length",
                default=23,
                description="K-mer length must be smaller or equal to the length of adapters.",
            )
            min_k = IntegerField(
                label="Minimum k-mer length at right end of reads used for trimming",
                default=11,
                disabled="preprocessing.adapters.length === 0 && preprocessing.custom_adapter_sequences.length === 0",
            )
            hamming_distance = IntegerField(
                label="Maximum Hamming distance for k-mers",
                default=1,
            )
            maxns = IntegerField(
                label="Max Ns after trimming",
                default=-1,
                description="If non-negative, reads with more Ns than this (after trimming) will be discarded.",
            )
            trim_quality = IntegerField(
                label="Quality below which to trim reads from the right end",
                default=10,
                description="Phred algorithm is used, which is more accurate than naive trimming.",
            )
            min_length = IntegerField(
                label="Minimum read length",
                default=20,
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

        class Quantification:
            """Quantification (Salmon)."""

            seq_bias = BooleanField(
                label="Perform sequence-specific bias correction",
                default=True,
                description="Perform sequence-specific bias correction.",
            )
            gc_bias = BooleanField(
                label="Perform fragment GC bias correction",
                required=False,
                description="Perform fragment GC bias correction. If single-end reads are selected "
                "as input in this workflow, it is recommended that you set this option to False. If you "
                "selected paired-end reads as input in this workflow, it is recommended that you set "
                "this option to True.",
            )
            consensus_slack = FloatField(
                label="Consensus slack",
                required=False,
                description="The amount of slack allowed in the quasi-mapping consensus mechanism. "
                "Normally, a transcript must cover all hits to be considered for mapping. "
                "If this is set to a fraction, X, greater than 0 (and in [0,1)), then a transcript "
                "can fail to cover up to (100 * X)% of the hits before it is discounted as a "
                "mapping candidate. The default value of this option is 0.2 in selective alignment mode "
                "and 0 otherwise.",
            )
            min_score_fraction = FloatField(
                label="Minimum alignment score fraction",
                default=0.65,
                description="The fraction of the optimal possible alignment score that a mapping "
                "must achieve in order to be considered valid - should be in (0,1].",
            )
            range_factorization_bins = IntegerField(
                label="Range factorization bins",
                default=4,
                description="Factorizes the likelihood used in quantification by adopting a "
                "new notion of equivalence classes based on the conditional probabilities with which "
                "fragments are generated from different transcripts. This is a more fine-grained "
                "factorization than the normal rich equivalence classes. The default value (4) "
                "corresponds to the default used in Zakeri et al. 2017 and larger values imply a more "
                "fine-grained factorization. If range factorization is enabled, a common value to "
                "select for this parameter is 4. A value of 0 signifies the use of basic rich "
                "equivalence classes.",
            )
            min_assigned_frag = IntegerField(
                label="Minimum number of assigned fragments",
                default=10,
                description="The minimum number of fragments that must be assigned to the "
                "transcriptome for quantification to proceed.",
            )
            num_bootstraps = IntegerField(
                label="--numBootstraps",
                description="Salmon has the ability to optionally "
                "compute bootstrapped abundance estimates. This is "
                "done by resampling (with replacement) from the counts "
                "assigned to the fragment equivalence classes, and then "
                "re-running the optimization procedure, either the EM or VBEM, "
                "for each such sample. The values of these different bootstraps "
                "allows us to assess technical variance in the main abundance "
                "estimates we produce. Such estimates can be useful for downstream "
                "(e.g. differential expression) tools that can make use of such "
                "uncertainty estimates. This option takes a positive integer that "
                "dictates the number of bootstrap samples to compute. The more samples "
                "computed, the better the estimates of varaiance, but the more "
                "computation (and time) required.",
                disabled="quantification.num_gibbs_samples",
                required=False,
            )
            num_gibbs_samples = IntegerField(
                label="--numGibbsSamples",
                description="Just as with the bootstrap procedure above, this option "
                "produces samples that allow us to estimate the variance in abundance "
                "estimates. However, in this case the samples are generated using posterior "
                "Gibbs sampling over the fragment equivalence classes rather than "
                "bootstrapping. We are currently analyzing these different approaches to "
                "assess the potential trade-offs in time / accuracy. The --numBootstraps "
                "and --numGibbsSamples options are mutually exclusive (i.e. in a given run, "
                "you must set at most one of these options to a positive integer.)",
                disabled="quantification.num_bootstraps",
                required=False,
            )

        class Downsampling:
            """Downsampling (Seqtk)."""

            n_reads = IntegerField(
                label="Number of reads",
                default=10000000,
                description="Number of reads to include in subsampling.",
            )
            seed = IntegerField(
                label="Number of reads",
                default=11,
                description="Using the same random seed makes reads subsampling reproducible "
                "in different environments.",
            )
            fraction = FloatField(
                label="Fraction of reads",
                required=False,
                range=[0.0, 1.0],
                description="Use the fraction of reads [0.0 - 1.0] from the orignal input file instead "
                "of the absolute number of reads. If set, this will override the 'Number of reads' "
                "input parameter.",
            )
            two_pass = BooleanField(
                label="2-pass mode",
                default=False,
                description="Enable two-pass mode when down-sampling. Two-pass mode is twice "
                "as slow but with much reduced memory usage.",
            )

        preprocessing = GroupField(
            Preprocessing,
            label="Preprocessing with BBDuk",
        )
        quantification = GroupField(
            Quantification,
            label="Quantification (Salmon)",
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
            bbduk_slug = "bbduk-single"
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            input_bbduk["operations"]["trim_pairs_evenly"] = True
            input_bbduk["operations"]["trim_by_overlap"] = True
            bbduk_slug = "bbduk-paired"
        else:
            self.error("Wrong reads input type.")

        preprocessing = Data.create(
            process=BioProcess.get_latest(bbduk_slug),
            input=input_bbduk,
            name=f"Trimmed ({inputs.reads.name})",
        )

        input_salmon = {
            "reads": preprocessing,
            "salmon_index": inputs.salmon_index,
            "annotation": inputs.annotation,
            "options": {
                "seq_bias": inputs.quantification.seq_bias,
                "min_score_fraction": inputs.quantification.min_score_fraction,
                "range_factorization_bins": inputs.quantification.range_factorization_bins,
                "min_assigned_frag": inputs.quantification.min_assigned_frag,
            },
        }

        if inputs.quantification.consensus_slack:
            input_salmon["options"][
                "consensus_slack"
            ] = inputs.quantification.consensus_slack

        if inputs.quantification.gc_bias:
            input_salmon["options"]["gc_bias"] = inputs.quantification.gc_bias

        if inputs.quantification.num_bootstraps:
            input_salmon["options"][
                "num_bootstraps"
            ] = inputs.quantification.num_bootstraps

        if inputs.quantification.num_gibbs_samples:
            input_salmon["options"][
                "num_gibbs_samples"
            ] = inputs.quantification.num_gibbs_samples

        quantification = Data.create(
            process=BioProcess.get_latest(slug="salmon-quant"),
            input=input_salmon,
            name=f"Quantified ({inputs.reads.name})",
        )

        input_seqtk = {
            "reads": preprocessing,
            "n_reads": inputs.downsampling.n_reads,
            "advanced": {
                "seed": inputs.downsampling.seed,
                "fraction": inputs.downsampling.fraction,
                "two_pass": inputs.downsampling.two_pass,
            },
        }

        if inputs.reads.type.startswith("data:reads:fastq:single:"):
            seqtk_slug = "seqtk-sample-single"
        elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
            seqtk_slug = "seqtk-sample-paired"
        else:
            self.error("Wrong reads input type.")

        downsampling = Data.create(
            process=BioProcess.get_latest(slug=seqtk_slug),
            input=input_seqtk,
            name=f"Subsampled ({inputs.reads.name})",
        )

        alignment_qc = Data.create(
            process=BioProcess.get_latest(slug="alignment-star"),
            input={
                "reads": downsampling,
                "genome": inputs.genome,
            },
            name=f"Aligned subset ({inputs.reads.name})",
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
                "alignment": alignment_qc,
            },
            name=f"Alignment summary ({inputs.reads.name})",
        )

        input_qorts = {
            "alignment": alignment_qc,
            "annotation": inputs.annotation,
            "options": {
                "stranded": "auto",
                "cdna_index": inputs.salmon_index,
                "n_reads": 5000000,
            },
        }

        qorts = Data.create(
            process=BioProcess.get_latest(slug="qorts-qc"),
            input=input_qorts,
            name=f"QoRTs QC report ({inputs.reads.name})",
        )

        input_multiqc = {
            "data": [
                inputs.reads,
                preprocessing,
                quantification,
                downsampling,
                alignment_qc,
                alignment_qc_rrna,
                alignment_qc_globin,
                idxstats,
                qorts,
            ]
        }

        Data.create(process=BioProcess.get_latest(slug="multiqc"), input=input_multiqc)

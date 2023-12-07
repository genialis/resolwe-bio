"""Beta Cut & Run workflow."""

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


class WorkflowCUTnRUN(Process):
    """Beta Cut & Run workflow.

    Analysis of samples processed for high resolution mapping of DNA binding sites using
    targeted nuclease strategy. The process is named CUT&RUN which stands for
    Cleavage Under Target and Release Using Nuclease. Workflow includes steps of
    trimming the reads with trimgalore, aligning them using bowtie2 to target species
    genome as well as a spike-in genome (optional). Aligned reads are processed to produce
    bigwig files to be viewed in a genome browser.

    """

    slug = "workflow-cutnrun-beta"
    name = "Beta Cut & Run workflow"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    version = "1.0.0"
    entity = {
        "type": "sample",
    }
    process_type = "data:workflow:cutnrun"
    category = "Pipeline"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq:paired",
            label="Input Reads (FASTQ)",
            description="Paired-end reads in FASTQ file.",
        )

        class TrimmingOptions:
            """Trimming options."""

            quality = IntegerField(
                label="Quality cutoff",
                description="Trim low-quality ends from reads based on phred score. "
                "Default: 20.",
                default=20,
            )
            nextseq = IntegerField(
                label="NextSeq/NovaSeq trim cutoff",
                description="NextSeq/NovaSeq-specific quality "
                "trimming. Trims also dark cycles appearing as "
                "high-quality G bases. This will set a specific "
                "quality cutoff, but qualities of G bases are ignored. "
                "This can not be used with Quality cutoff and will "
                "override it.",
                required=False,
            )
            min_length = IntegerField(
                label="Minimum length after trimming",
                description="Discard reads that became shorter than "
                "selected length because of either quality or adapter "
                "trimming. Both reads of a read-pair need to be longer "
                "than specified length to be printed out to validated "
                "paired-end files. If only one read became too short "
                "there is the possibility of keeping such unpaired "
                "single-end reads with Retain unpaired. A value of 0 "
                "disables filtering based on length. Default: 20.",
                default=20,
            )

            class AdapterOptions:
                """Adapter trimming options."""

                adapter_1 = ListField(
                    StringField(),
                    label="Read 1 adapter sequence",
                    description="Adapter sequences to be trimmed. "
                    "Also see universal adapters field for predefined "
                    "adapters. This is mutually exclusive with read 1 "
                    "adapters file and universal adapters.",
                    required=False,
                    default=[],
                )
                adapter_2 = ListField(
                    StringField(),
                    label="Read 2 adapter sequence",
                    description="Optional adapter sequence to be trimmed "
                    "off read 2 of paired-end files. This is mutually "
                    "exclusive with read 2 adapters file and universal "
                    "adapters.",
                    required=False,
                    default=[],
                )
                adapter_file_1 = DataField(
                    "seq:nucleotide",
                    label="Read 1 adapters file",
                    description="This is mutually exclusive with read 1 "
                    "adapters and universal adapters.",
                    required=False,
                )
                adapter_file_2 = DataField(
                    "seq:nucleotide",
                    label="Read 2 adapters file",
                    description="This is mutually exclusive with read 2 "
                    "adapters and universal adapters.",
                    required=False,
                )
                universal_adapter = StringField(
                    label="Universal adapters",
                    description="Instead of default detection use specific "
                    "adapters. Use 13bp of the Illumina universal adapter, "
                    "12bp of the Nextera adapter or 12bp of the Illumina "
                    "Small RNA 3' Adapter. Selecting to trim smallRNA "
                    "adapters will also lower the min length value to 18bp. "
                    "If the smallRNA libraries are paired-end then read 2 "
                    "adapter will be set to the Illumina small RNA 5' "
                    "adapter automatically (GATCGTCGGACT) unless defined "
                    "explicitly. This is mutually exclusive with manually "
                    "defined adapters and adapter files.",
                    choices=[
                        ("--illumina", "Illumina"),
                        ("--nextera", "Nextera"),
                        ("--small_rna", "Illumina small RNA"),
                    ],
                    required=False,
                )
                stringency = IntegerField(
                    label="Overlap with adapter sequence required to trim ",
                    description="Defaults to a very stringent setting of "
                    "1, i.e. even a single base pair of overlapping "
                    "sequence will be trimmed of the 3' end of any read. "
                    "Default: 1.",
                    default=1,
                )
                error_rate = FloatField(
                    label="Maximum allowed error rate",
                    description="Number of errors divided by the length of "
                    "the matching region. Default: 0.1.",
                    default=0.1,
                )

            adapter_options = GroupField(
                AdapterOptions, label="Adapter trimming options"
            )

        class AlignmentOptions:
            """Alignment options."""

            genome = DataField(
                data_type="index:bowtie2",
                label="Species genome",
            )
            spike_in_genome = DataField(
                data_type="index:bowtie2",
                label="Spike-in genome",
                required=False,
                disabled="normalization_options.skip_norm == True",
            )
            alignment_mode = StringField(
                label="Alignment mode",
                choices=[
                    ("--local", "Local"),
                    ("--end-to-end", "End-to-end"),
                ],
                default="--end-to-end",
                description="Local: Some characters may be omitted ('soft clipped') from the "
                "ends in order to achieve the greatest possible alignment score. "
                "End-to-end: Option without any trimming (or 'soft clipping') of bases from either end. "
                "This option is enabled by default and is suitable if reads have been clipped beforehand.",
            )
            speed = StringField(
                label="Speed vs. Sensitivity",
                choices=[
                    ("--very-fast", "Very fast"),
                    ("--fast", "Fast"),
                    ("--sensitive", "Sensitive"),
                    ("--very-sensitive", "Very sensitive"),
                ],
                default="--very-sensitive",
                description="A quick setting for aligning fast or accurately."
                "This option is a shortcut for parameters as follows:"
                "For End-to-end:"
                "Very fast                   -D 5 -R 1 -N 0 -L 22 -i S,0,2.50"
                "Fast                        -D 10 -R 2 -N 0 -L 22 -i S,0,2.50"
                "Sensitive                   -D 15 -R 2 -N 0 -L 22 -i S,1,1.15"
                "Very sensitive (default)    -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"
                "For Local:"
                "Very fast                   -D 5 -R 1 -N 0 -L 25 -i S,1,2.00"
                "Fast                        -D 10 -R 2 -N 0 -L 22 -i S,1,1.75"
                "Sensitive                   -D 15 -R 2 -N 0 -L 20 -i S,1,0.75"
                "Very sensitive              -D 20 -R 3 -N 0 -L 20 -i S,1,0.50",
            )

            class PEOptions:
                """Paired-end options."""

                dovetail = BooleanField(
                    label="Dovetail",
                    default=True,
                    description="If the mates “dovetail”, that is if one mate "
                    "alignment extends past the beginning of the other such that "
                    "the wrong mate begins upstream. Consider "
                    "that to be concordant. Default: True.",
                )
                rep_se = BooleanField(
                    label="Report single ended",
                    default=False,
                    description="If paired alignment can not be found Bowtie2 tries "
                    "to find alignments for the individual mates. Default: False.",
                )
                minins = IntegerField(
                    label="Minimal distance",
                    default=10,
                    description="The minimum fragment length (--minins) for valid paired-end alignments. "
                    "Default: 10.",
                )
                maxins = IntegerField(
                    label="Maximal distance",
                    default=700,
                    description="The maximum fragment length (--maxins) for valid paired-end alignments. "
                    "Default: 700.",
                )
                discordantly = BooleanField(
                    label="Report discordantly matched read",
                    default=False,
                    description="If both mates have unique alignments, but the alignments do "
                    "not match paired-end expectations (orientation and relative distance) then "
                    "alignment will be reported. Useful for detecting structural variations. "
                    "Default: False.",
                )

            class OutputOptions:
                """Output options."""

                no_unal = BooleanField(
                    label="Suppress SAM records for unaligned reads",
                    default=True,
                    description="When enabled, suppress SAM records for unaligned reads. Default: True.",
                )

            pe_options = GroupField(PEOptions, label="Paired-end options")
            output_options = GroupField(OutputOptions, label="Output options")

        class NormalizationOptions:
            """Spike-in normalization options."""

            skip_norm = BooleanField(
                label="Skip normalization",
                default=False,
                description="Skip the spike-in normalization step of BigWig output. "
                "Use this if you don't provide a spike-in. Default: False.",
            )
            scale = FloatField(
                label="Scale factor",
                default=10000,
                description="Magnitude of the scale factor. The scaling factor is calculated by "
                "dividing the scale with the number of features in BEDPE (scale/(number of features)). "
                "Default: 10000.",
                disabled="normalization_options.skip_norm == True",
            )

        class SubsampleOptions:
            """Subsampling options."""

            downsample_reads = BooleanField(
                label="Downsample reads",
                default=True,
                description="Option to downsample reads before trimming. Default: True",
            )

            n_reads = IntegerField(
                label="Number of reads to subsample",
                default=10000000,
                description="Number of reads to subsample from the input FASTQ file. Default: 10M.",
                disabled="subsample_options.downsample_reads == False",
            )

        class DeduplicationOptions:
            """Deduplication options."""

            remove_duplicates = BooleanField(
                label="Remove duplicates",
                default=False,
                description="Option on how to handle duplicate reads. "
                "True: Mark and remove duplicate reads. "
                "False: Only mark duplicate reads. "
                "Note that this option is only available for species genome. "
                "In case of spike-in genome, duplicate reads are always removed. "
                "Default: False.",
            )

        trimming_options = GroupField(TrimmingOptions, label="Trimming options")
        alignment_options = GroupField(AlignmentOptions, label="Alignment options")
        normalization_options = GroupField(
            NormalizationOptions, label="Spike-in normalization options"
        )
        subsample_options = GroupField(SubsampleOptions, label="Subsample options")
        deduplication_options = GroupField(
            DeduplicationOptions, label="Deduplication options"
        )

    class Output:
        """Output fields."""

        # Workflows do not have output fields.

    def run(self, inputs, outputs):
        """Run the workflow."""

        if inputs.subsample_options.downsample_reads:
            downsampled_reads = Data.create(
                process=BioProcess.get_latest(slug="seqtk-sample-paired"),
                input={
                    "reads": inputs.reads,
                    "n_reads": inputs.subsample_options.n_reads,
                },
                name=f"Subset ({inputs.reads.name})",
            )
            reads = downsampled_reads
        else:
            reads = inputs.reads

        input_trimgalore = {
            "reads": reads,
            "quality_trim": {
                "quality": inputs.trimming_options.quality,
                "min_length": inputs.trimming_options.min_length,
            },
            "adapter_trim": {
                "stringency": inputs.trimming_options.adapter_options.stringency,
                "error_rate": inputs.trimming_options.adapter_options.error_rate,
            },
        }

        if inputs.trimming_options.nextseq:
            input_trimgalore["quality_trim"][
                "nextseq"
            ] = inputs.trimming_options.nextseq
        if inputs.trimming_options.adapter_options.adapter_1:
            input_trimgalore["adapter_trim"][
                "adapter"
            ] = inputs.trimming_options.adapter_options.adapter_1
        if inputs.trimming_options.adapter_options.adapter_2:
            input_trimgalore["adapter_trim"][
                "adapter_2"
            ] = inputs.trimming_options.adapter_options.adapter_2
        if inputs.trimming_options.adapter_options.adapter_file_1:
            input_trimgalore["adapter_trim"][
                "adapter_file_1"
            ] = inputs.trimming_options.adapter_options.adapter_file_1
        if inputs.trimming_options.adapter_options.adapter_file_2:
            input_trimgalore["adapter_trim"][
                "adapter_file_2"
            ] = inputs.trimming_options.adapter_options.adapter_file_2
        if inputs.trimming_options.adapter_options.universal_adapter:
            input_trimgalore["adapter_trim"][
                "universal_adapter"
            ] = inputs.trimming_options.adapter_options.universal_adapter

        preprocessed_reads = Data.create(
            process=BioProcess.get_latest(slug="trimgalore-paired"),
            input=input_trimgalore,
            name=f"Trimmed ({inputs.reads.name})",
        )

        input_alignment = {
            "reads": preprocessed_reads,
            "genome": inputs.alignment_options.genome,
            "mode": inputs.alignment_options.alignment_mode,
            "speed": inputs.alignment_options.speed,
            "PE_options": {
                "dovetail": inputs.alignment_options.pe_options.dovetail,
                "rep_se": inputs.alignment_options.pe_options.rep_se,
                "minins": inputs.alignment_options.pe_options.minins,
                "maxins": inputs.alignment_options.pe_options.maxins,
                "discordantly": inputs.alignment_options.pe_options.discordantly,
            },
            "output_opts": {
                "no_unal": inputs.alignment_options.output_options.no_unal,
            },
        }

        species_alignment = Data.create(
            process=BioProcess.get_latest(slug="alignment-bowtie2"),
            input=input_alignment,
            name=f"Alignment species ({inputs.reads.name})",
        )

        input_deduplicate = {
            "bam": species_alignment,
            "remove_duplicates": True
            if inputs.deduplication_options.remove_duplicates == True
            else False,
        }

        deduplicate_alignment = Data.create(
            process=BioProcess.get_latest(
                slug="markduplicates",
            ),
            input=input_deduplicate,
            name=f"Deduplicated alignment species ({inputs.reads.name})",
        )

        if not inputs.normalization_options.skip_norm:
            input_alignment_spikein = {
                "genome": inputs.alignment_options.spike_in_genome,
                "reads": reads,
                "mode": inputs.alignment_options.alignment_mode,
                "speed": inputs.alignment_options.speed,
                "PE_options": {
                    "dovetail": inputs.alignment_options.pe_options.dovetail,
                    "rep_se": inputs.alignment_options.pe_options.rep_se,
                    "minins": inputs.alignment_options.pe_options.minins,
                    "maxins": inputs.alignment_options.pe_options.maxins,
                    "discordantly": inputs.alignment_options.pe_options.discordantly,
                },
                "output_opts": {
                    "no_unal": inputs.alignment_options.output_options.no_unal,
                },
            }

            spikein_alignment = Data.create(
                process=BioProcess.get_latest(slug="alignment-bowtie2"),
                input=input_alignment_spikein,
                name=f"Alignment spike-in ({inputs.reads.name})",
            )

            deduplicate_alignment_spikein = Data.create(
                process=BioProcess.get_latest(
                    slug="markduplicates",
                ),
                input={"bam": spikein_alignment, "remove_duplicates": True},
                name=f"Deduplicated alignment spike-in ({inputs.reads.name})",
            )

            spikein_normfactor = Data.create(
                process=BioProcess.get_latest(
                    slug="bedtools-bamtobed",
                ),
                input={"alignment": deduplicate_alignment_spikein},
            )

            input_normalization = {
                "alignment": deduplicate_alignment,
                "bedpe": spikein_normfactor,
                "scale": inputs.normalization_options.scale,
            }

        else:
            input_normalization = {
                "alignment": deduplicate_alignment,
            }

        Data.create(
            process=BioProcess.get_latest(
                slug="calculate-bigwig",
            ),
            input=input_normalization,
        )

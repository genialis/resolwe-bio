"""ATAC-Seq workflow."""

from resolwe.process import (
    BooleanField,
    Data,
    DataField,
    FloatField,
    GroupField,
    IntegerField,
    Process,
    StringField,
)
from resolwe.process.models import Process as BioProcess


class WorkflowATACSeq(Process):
    """Beta ATAC-Seq workflow.

    This ATAC-seq pipeline follows the official ENCODE DCC pipeline with additional modifications.

    First, reads are optionally downsampled and then aligned to a genome using
    [Bowtie2](http://bowtie-bio.sourceforge.net/index.shtml) aligner. Next, pre-peakcall QC
    metrics are calculated. QC report contains ENCODE 3 proposed QC metrics --
    [NRF](https://www.encodeproject.org/data-standards/terms/),
    [PBC bottlenecking coefficients, NSC, and RSC](https://genome.ucsc.edu/ENCODE/qualityMetrics.html#chipSeq).
    The peaks are called using [MACS2](https://github.com/taoliu/MACS/).
    The post-peakcall QC report includes additional QC metrics -- number of peaks,
    fraction of reads in peaks (FRiP), number of reads in peaks, and if promoter
    regions BED file is provided, number of reads in promoter regions, fraction of
    reads in promoter regions, number of peaks in promoter regions, and fraction of reads in promoter regions.
    Finally, a bigWig file is produced along with the ChipQC and MultiQC reports.
    """

    slug = "workflow-atac-seq-beta"
    name = "ATAC-Seq (Beta)"
    requirements = {
        "expression-engine": "jinja",
    }
    data_name = "{{ reads|name|default('?') }}"
    version = "0.1.0"
    entity = {"type": "sample"}
    process_type = "data:workflow:atacseq"
    category = "ATAC-seq"

    class Input:
        """Input fields."""

        reads = DataField(
            data_type="reads:fastq",
            label="Select sample(s)",
        )
        genome = DataField(
            data_type="index:bowtie2",
            label="Genome",
        )
        promoter = DataField(
            data_type="bed",
            label="Promoter regions BED file",
            required=False,
            description="BED file containing promoter regions (TSS+-1000 bp for example). "
            "Needed to get the number of peaks and reads mapped to promoter regions.",
        )

        class DownsamplingOptions:
            """Downsampling options."""

            downsample_reads = BooleanField(
                label="Downsample reads",
                default=False,
                description="If selected, reads will be downsampled to the specified number of reads "
                "before alignment.",
            )
            n_reads = IntegerField(
                label="Number of reads",
                default=10000000,
                disabled="!downsampling_options.downsample_reads",
                description="Number of reads to include in downsampling.",
            )
            two_pass = BooleanField(
                label="2-pass mode",
                default=False,
                disabled="!downsampling_options.downsample_reads",
                description="Enable two-pass mode when down-sampling. "
                "Two-pass mode is twice as slow but with much reduced memory.",
            )

        downsampling_options = GroupField(DownsamplingOptions, label="Downsampling")

        class AlignmentOptions:
            """Alignment (Bowtie2) options."""

            mode = StringField(
                label="Alignment mode",
                default="--local",
                choices=[
                    ("--end-to-end", "end to end mode"),
                    ("--local", "local"),
                ],
                description="End to end: Bowtie 2 requires that the entire read align from one end to the other, "
                "without any trimming (or 'soft clipping') of characters from either end. "
                "Local: Bowtie 2 does not require that the entire read align from one end to the other. "
                "Rather, some characters may be omitted ('soft clipped') from the ends in order to "
                "achieve the greatest possible alignment score.",
            )
            speed = StringField(
                label="Speed vs. Sensitivity",
                default="--sensitive",
                choices=[
                    ("--very-fast", "Very fast"),
                    ("--fast", "Fast"),
                    ("--sensitive", "Sensitive"),
                    ("--very-sensitive", "Very sensitive"),
                ],
            )

            class PEOptions:
                """Paired end alignment options."""

                use_se = BooleanField(
                    label="Map as single-ended (for paired-end reads only)",
                    default=False,
                    description="If this option is selected paired-end reads will be mapped as single-ended and "
                    "other paired-end options are ignored.",
                )
                discordantly = BooleanField(
                    label="Report discordantly matched read",
                    default=True,
                    description="If both mates have unique alignments, but the alignments do not match paired-end "
                    "expectations (orientation and relative distance) then alignment will be reported. "
                    "Useful for detecting structural variations.",
                )
                rep_se = BooleanField(
                    label="Report single ended",
                    default=True,
                    description="If paired alignment can not be found Bowtie2 tries to find alignments for the "
                    "individual mates.",
                )
                minins = IntegerField(
                    label="Minimal distance",
                    default=0,
                    description="The minimum fragment length for valid paired-end alignments. "
                    "0 imposes no minimum.",
                )
                maxins = IntegerField(
                    label="Maximal distance",
                    default=2000,
                    description="The maximum fragment length for valid paired-end alignments.",
                )

            pe_options = GroupField(
                PEOptions,
                label="Paired end alignment options",
            )

            class TrimmingOptions:
                """Trimming options."""

                trim_5 = IntegerField(
                    label="Bases to trim from 5'",
                    default=0,
                    description="Number of bases to trim from the 5' (left) end of each read before alignment.",
                )
                trim_3 = IntegerField(
                    label="Bases to trim from 3'",
                    default=0,
                    description="Number of bases to trim from the 3' (right) end of each read before alignment.",
                )
                trim_iter = IntegerField(
                    label="Iterations",
                    default=0,
                    description="Number of iterations.",
                )
                trim_nucl = IntegerField(
                    label="Bases to trim",
                    default=2,
                    description="Number of bases to trim from 3' end in each iteration.",
                )

            trimming_options = GroupField(TrimmingOptions, label="Iterative trimming")

            class ReportingOptions:
                """Reporting options."""

                rep_mode = StringField(
                    label="Report mode",
                    default="def",
                    choices=[
                        ("def", "Default mode"),
                        ("k", "-k mode"),
                        ("a", "-a mode (very slow)"),
                    ],
                    description="Default mode: search for multiple alignments, report the best one; "
                    "-k mode: search for one or more alignments, report each; "
                    "-a mode: search for and report all alignments",
                )
                k_reports = IntegerField(
                    label="Number of reports (for -k mode only)",
                    default=5,
                    description="Searches for at most X distinct, valid alignments for each read. The search "
                    "terminates when it can't find more distinct valid alignments, or when it finds X, "
                    "whichever happens first.",
                )

            reporting_options = GroupField(ReportingOptions, label="Reporting")

        alignment_options = GroupField(AlignmentOptions, label="Alignment (Bowtie2)")

        class PrepeakQCOptions:
            """Pre-peak QC settings."""

            q_threshold = IntegerField(
                label="Quality filtering threshold",
                default=30,
            )
            n_sub = IntegerField(
                label="Number of reads to subsample",
                default=25000000,
            )
            tn5 = BooleanField(
                label="Tn5 shifting",
                default=True,
                description="Tn5 transposon shifting. Shift reads on '+' strand by 4 bp "
                "and reads on '-' strand by 5 bp.",
            )
            shift = IntegerField(
                label="User-defined cross-correlation peak strandshift",
                default=0,
                description="If defined, SPP tool will not try to estimate fragment length "
                "but will use the given value as fragment length.",
            )

        prepeakqc_options = GroupField(PrepeakQCOptions, label="Pre-peak QC settings")

        class Macs2Options:
            """MACS2 options."""

            tagalign = BooleanField(
                label="Use tagAlign files",
                default=True,
                description="Use filtered tagAlign files as case (treatment) and control "
                "(background) samples. If extsize parameter is not set, run MACS "
                "using input's estimated fragment length.",
            )
            duplicates = StringField(
                label="Number of duplicates",
                required=False,
                hidden="macs2_options.tagalign",
                choices=[
                    ("1", "1"),
                    ("auto", "auto"),
                    ("all", "all"),
                ],
                description="It controls the MACS behavior towards duplicate tags at the exact same location -- the "
                "same coordination and the same strand. The 'auto' option makes MACS calculate the "
                "maximum tags at the exact same location based on binomal distribution using 1e-5 as "
                "pvalue cutoff and the 'all' option keeps all the tags. If an integer is given, at most "
                "this number of tags will be kept at the same location. The default is to keep one tag "
                "at the same location.",
            )
            duplicates_prepeak = StringField(
                label="Number of duplicates",
                required=False,
                default="all",
                hidden="!macs2_options.tagalign",
                choices=[
                    ("1", "1"),
                    ("auto", "auto"),
                    ("all", "all"),
                ],
                description="It controls the MACS behavior towards duplicate tags at the exact same location -- the "
                "same coordination and the same strand. The 'auto' option makes MACS calculate the "
                "maximum tags at the exact same location based on binomal distribution using 1e-5 as "
                "pvalue cutoff and the 'all' option keeps all the tags. If an integer is given, at most "
                "this number of tags will be kept at the same location. The default is to keep one tag "
                "at the same location.",
            )
            qvalue = FloatField(
                label="Q-value cutoff",
                required=False,
                disabled="macs2_options.pvalue && macs2_options.pvalue_prepeak",
                description="The q-value (minimum FDR) cutoff to call significant regions. Q-values "
                "are calculated from p-values using Benjamini-Hochberg procedure.",
            )
            pvalue = FloatField(
                label="P-value cutoff",
                required=False,
                hidden="macs2_options.tagalign",
                disabled="macs2_options.qvalue",
                description="The p-value cutoff. If specified, MACS2 will use p-value instead of q-value cutoff.",
            )
            pvalue_prepeak = FloatField(
                label="P-value cutoff",
                default=0.01,
                hidden="!macs2_options.tagalign || macs2_options.qvalue",
                disabled="macs2_options.qvalue",
                description="The p-value cutoff. If specified, MACS2 will use p-value instead of q-value cutoff.",
            )
            cap_num = IntegerField(
                label="Cap number of peaks by taking top N peaks",
                default=300000,
                disabled="macs2_options.broad",
                description="To keep all peaks set value to 0.",
            )
            mfold_lower = IntegerField(
                label="MFOLD range (lower limit)",
                required=False,
                description="This parameter is used to select the regions within MFOLD range of high-confidence "
                "enrichment ratio against background to build model. The regions must be lower than "
                "upper limit, and higher than the lower limit of fold enrichment. DEFAULT:10,30 means "
                "using all regions not too low (>10) and not too high (<30) to build paired-peaks "
                "model. If MACS can not find more than 100 regions to build model, it will use the "
                "--extsize parameter to continue the peak detection ONLY if --fix-bimodal is set.",
            )
            mfold_upper = IntegerField(
                label="MFOLD range (upper limit)",
                required=False,
                description="This parameter is used to select the regions within MFOLD range of high-confidence "
                "enrichment ratio against background to build model. The regions must be lower than "
                "upper limit, and higher than the lower limit of fold enrichment. DEFAULT:10,30 means "
                "using all regions not too low (>10) and not too high (<30) to build paired-peaks "
                "model. If MACS can not find more than 100 regions to build model, it will use the "
                "--extsize parameter to continue the peak detection ONLY if --fix-bimodal is set.",
            )
            slocal = IntegerField(
                label="Small local region",
                required=False,
                description="Slocal and llocal parameters control which two levels of regions will be checked "
                "around the peak regions to calculate the maximum lambda as local lambda. By default, "
                "MACS considers 1000 bp for small local region (--slocal), and 10000 bp for large local "
                "region (--llocal) which captures the bias from a long range effect like an open "
                "chromatin domain. You can tweak these according to your project. Remember that if the "
                "region is set too small, a sharp spike in the input data may kill the significant peak.",
            )
            llocal = IntegerField(
                label="Large local region",
                required=False,
                description="Slocal and llocal parameters control which two levels of regions will be checked "
                "around the peak regions to calculate the maximum lambda as local lambda. By default, "
                "MACS considers 1000 bp for small local region (--slocal), and 10000 bp for large local "
                "region (--llocal) which captures the bias from a long range effect like an open "
                "chromatin domain. You can tweak these according to your project. Remember that if the "
                "region is set too small, a sharp spike in the input data may kill the significant peak.",
            )
            extsize = IntegerField(
                label="extsize",
                default=150,
                description="While '--nomodel' is set, MACS uses this parameter to extend reads in 5'->3' direction "
                "to fix-sized fragments. For example, if the size of binding region for your "
                "transcription factor is 200 bp, and you want to bypass the model building by MACS, "
                "this parameter can be set as 200. This option is only valid when --nomodel is set or "
                "when MACS fails to build model and --fix-bimodal is on.",
            )
            shift = IntegerField(
                label="Shift",
                default=-75,
                description="Note, this is NOT the legacy --shiftsize option which is replaced by --extsize! You "
                "can set an arbitrary shift in bp here. Please Use discretion while setting it other "
                "than default value (0). When --nomodel is set, MACS will use this value to move "
                "cutting ends (5') then apply --extsize from 5' to 3' direction to extend them to "
                "fragments. When this value is negative, ends will be moved toward 3'->5' direction,"
                "otherwise 5'->3' direction. Recommended to keep it as default 0 for ChIP-Seq datasets, "
                "or -1 * half of EXTSIZE together with --extsize option for detecting enriched cutting "
                "loci such as certain DNAseI-Seq datasets. Note, you can't set values other than 0 if "
                "format is BAMPE for paired-end data. Default is 0.",
            )
            band_width = IntegerField(
                label="Band width",
                required=False,
                description="The band width which is used to scan the genome ONLY for model building. You can set "
                "this parameter as the sonication fragment size expected from wet experiment. The "
                "previous side effect on the peak detection process has been removed. So this parameter "
                "only affects the model building.",
            )
            nolambda = BooleanField(
                label="Use background lambda as local lambda",
                default=False,
                description="With this flag on, MACS will use the background lambda as local lambda. This means "
                "MACS will not consider the local bias at peak candidate regions.",
            )
            fix_bimodal = BooleanField(
                label="Turn on the auto paired-peak model process",
                default=False,
                description="Turn on the auto paired-peak model process. If it's set, when MACS failed "
                "to build paired model, it will use the nomodel settings, the '--extsize' parameter "
                "to extend each tag. If set, MACS will be terminated if paired-peak model has failed.",
            )
            nomodel = BooleanField(
                label="Bypass building the shifting model",
                default=False,
                hidden="macs2_options.tagalign",
                description="While on, MACS will bypass building the shifting model.",
            )
            nomodel_prepeak = BooleanField(
                label="Bypass building the shifting model",
                default=True,
                hidden="!macs2_options.tagalign",
                description="While on, MACS will bypass building the shifting model.",
            )
            down_sample = BooleanField(
                label="Down-sample",
                default=False,
                description="When set to true, random sampling method will scale down the bigger sample. By default, MACS "
                "uses linear scaling. This option will make the results unstable and irreproducible "
                "since each time, random reads would be selected, especially the numbers (pileup, "
                "pvalue, qvalue) would change.",
            )
            bedgraph = BooleanField(
                label="Save fragment pileup and control lambda",
                default=True,
                description="If this flag is on, MACS will store the fragment pileup, control lambda, -log10pvalue "
                "and -log10qvalue scores in bedGraph files. The bedGraph files will be stored in "
                "current directory named NAME+'_treat_pileup.bdg' for treatment data, "
                "NAME+'_control_lambda.bdg' for local lambda values from control, "
                "NAME+'_treat_pvalue.bdg' for Poisson pvalue scores (in -log10(pvalue) form), and "
                "NAME+'_treat_qvalue.bdg' for q-value scores from Benjamini-Hochberg-Yekutieli "
                "procedure.",
            )
            spmr = BooleanField(
                label="Save signal per million reads for fragment pileup profiles",
                default=True,
                disabled="macs2_options.bedgraph === false",
            )
            call_summits = BooleanField(
                label="Call summits",
                default=True,
                description="MACS will now reanalyze the shape of signal profile (p or q-score depending on cutoff "
                "setting) to deconvolve subpeaks within each peak called from general procedure. It's "
                "highly recommended to detect adjacent binding events. While used, the output subpeaks "
                "of a big peak region will have the same peak boundaries, and different scores and peak "
                "summit positions.",
            )
            broad = BooleanField(
                label="Composite broad regions",
                default=False,
                disabled="macs2_options.call_summits === true",
                description="When this flag is on, MACS will try to composite broad regions in BED12 (a "
                "gene-model-like format) by putting nearby highly enriched regions into a broad region "
                "with loose cutoff. The broad region is controlled by another cutoff through "
                "--broad-cutoff. The maximum length of broad region length is 4 times of d from MACS.",
            )
            broad_cutoff = FloatField(
                label="Broad cutoff",
                required=False,
                disabled="macs2_options.call_summits === true || macs2_options.broad !== true",
                description="Cutoff for broad region. This option is not available unless --broad is set. If -p is "
                "set, this is a p-value cutoff, otherwise, it's a q-value cutoff. DEFAULT = 0.1",
            )

        macs2_options = GroupField(Macs2Options, label="MACS2 options")

        class BigWigOptions:
            """bigWig options."""

            bin_size = IntegerField(
                label="Bin size",
                default=50,
                description="Size of the bins (in bp) for the output bigWig file. "
                "A smaller bin size value will result in a higher resolution of "
                "the coverage track but also in a larger file size.",
            )

        bigwig_options = GroupField(BigWigOptions, label="bigWig options")

    class Output:
        """Output fields."""

    def run(self, inputs, outputs):
        """Run the workflow."""

        reads = inputs.reads

        if inputs.downsampling_options.downsample_reads:
            seqtk_input = {
                "reads": inputs.reads,
                "n_reads": inputs.downsampling_options.n_reads,
                "advanced": {
                    "two_pass": inputs.downsampling_options.two_pass,
                },
            }

            if inputs.reads.type.startswith("data:reads:fastq:single:"):
                seqtk_slug = "seqtk-sample-single"
            elif inputs.reads.type.startswith("data:reads:fastq:paired:"):
                seqtk_slug = "seqtk-sample-paired"

            reads = Data.create(
                process=BioProcess.get_latest(slug=seqtk_slug),
                input=seqtk_input,
                name=f"Subset ({inputs.reads.name})",
            )

        alignment_input = {
            "genome": inputs.genome,
            "reads": reads,
            "mode": inputs.alignment_options.mode,
            "speed": inputs.alignment_options.speed,
            "PE_options": {
                "use_se": inputs.alignment_options.pe_options.use_se,
                "discordantly": inputs.alignment_options.pe_options.discordantly,
                "rep_se": inputs.alignment_options.pe_options.rep_se,
                "minins": inputs.alignment_options.pe_options.minins,
                "maxins": inputs.alignment_options.pe_options.maxins,
            },
            "start_trimming": {
                "trim_5": inputs.alignment_options.trimming_options.trim_5,
                "trim_3": inputs.alignment_options.trimming_options.trim_3,
            },
            "trimming": {
                "trim_iter": inputs.alignment_options.trimming_options.trim_iter,
                "trim_nucl": inputs.alignment_options.trimming_options.trim_nucl,
            },
            "reporting": {
                "rep_mode": inputs.alignment_options.reporting_options.rep_mode,
                "k_reports": inputs.alignment_options.reporting_options.k_reports,
            },
        }

        alignment = Data.create(
            process=BioProcess.get_latest(slug="alignment-bowtie2"),
            input=alignment_input,
            name=f"Aligned ({inputs.reads.name})",
        )

        bigwig_input = {
            "alignment": alignment,
            "bin_size": inputs.bigwig_options.bin_size,
        }
        Data.create(
            process=BioProcess.get_latest(slug="calculate-bigwig"),
            input=bigwig_input,
            name=f"Coverage ({inputs.reads.name})",
        )

        macs_input = {
            "case": alignment,
            "tagalign": inputs.macs2_options.tagalign,
            "prepeakqc_settings": {
                "q_threshold": inputs.prepeakqc_options.q_threshold,
                "n_sub": inputs.prepeakqc_options.n_sub,
                "tn5": inputs.prepeakqc_options.tn5,
                "shift": inputs.prepeakqc_options.shift,
            },
            "settings": {
                "duplicates": inputs.macs2_options.duplicates,
                "duplicates_prepeak": inputs.macs2_options.duplicates_prepeak,
                "qvalue": inputs.macs2_options.qvalue,
                "pvalue": inputs.macs2_options.pvalue,
                "pvalue_prepeak": inputs.macs2_options.pvalue_prepeak,
                "cap_num": inputs.macs2_options.cap_num,
                "mfold_lower": inputs.macs2_options.mfold_lower,
                "mfold_upper": inputs.macs2_options.mfold_upper,
                "slocal": inputs.macs2_options.slocal,
                "llocal": inputs.macs2_options.llocal,
                "extsize": inputs.macs2_options.extsize,
                "shift": inputs.macs2_options.shift,
                "nolambda": inputs.macs2_options.nolambda,
                "fix_bimodal": inputs.macs2_options.fix_bimodal,
                "nomodel": inputs.macs2_options.nomodel,
                "nomodel_prepeak": inputs.macs2_options.nomodel_prepeak,
                "down_sample": inputs.macs2_options.down_sample,
                "bedgraph": inputs.macs2_options.bedgraph,
                "spmr": inputs.macs2_options.spmr,
                "call_summits": inputs.macs2_options.call_summits,
                "broad": inputs.macs2_options.broad,
                "broad_cutoff": inputs.macs2_options.broad_cutoff,
            },
        }

        if inputs.promoter:
            macs_input["promoter"] = inputs.promoter

        macs2 = Data.create(
            process=BioProcess.get_latest(slug="macs2-callpeak"),
            input=macs_input,
            name=f"MACS2 peaks ({inputs.reads.name})",
        )

        input_chipqc = {
            "alignment": alignment,
            "peaks": macs2,
        }

        chipqc = Data.create(
            process=BioProcess.get_latest(
                slug="chipqc",
            ),
            input=input_chipqc,
        )

        Data.create(
            process=BioProcess.get_latest(slug="multiqc"),
            input={"data": [reads, alignment, macs2, chipqc]},
            name=f"MultiQC ({inputs.reads.name})",
        )

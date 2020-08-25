"""Calculate quality metrics for ChIPseq data.."""

from pathlib import Path
from shutil import copyfile

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DirField,
    FileField,
    GroupField,
    IntegerField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class ChipQC(Process):
    """Calculate quality control metrics for ChIP-seq samples.

    The analysis is based on ChIPQC package which computs a variety of
    quality control metrics and statistics, and provides plots and
    a report for assessment of experimental data for further analysis.
    """

    slug = "chipqc"
    name = "ChipQC"
    process_type = "data:chipqc"
    version = "1.0.3"
    category = "ChIP-Seq:QC report"
    data_name = '{{ alignment|sample_name|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/chipseq:4.1.3"}},
        "resources": {
            "cores": 8,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process ChipQC."""

        alignment = DataField(
            data_type="alignment:bam",
            label="Aligned reads",
        )
        peaks = DataField(
            data_type="chipseq:callpeak",
            label="Called peaks",
        )
        blacklist = DataField(
            data_type="bed",
            label="Blacklist regions",
            description="BED file containing genomic regions that should be "
            "excluded from the analysis.",
            required=False,
        )
        calculate_enrichment = BooleanField(
            label="Calculate enrichment",
            description="Calculate enrichment of signal in known genomic "
            "annotation. By default annotation is provided from "
            "the TranscriptDB package specified by genome bulid "
            "which should match one of the supported annotations "
            "(hg19, hg38, hg18, mm10, mm9, rn4, ce6, dm3). If "
            "annotation is not supported the analysis is skipped.",
            default=False,
        )

        class Advanced:
            """Add advanced list of options."""

            quality_threshold = IntegerField(
                label="Mapping quality threshold",
                description="Only reads with mapping quality scores above "
                "this threshold will be used for some statistics.",
                default=15,
            )
            profile_window = IntegerField(
                label="Window size",
                description="An integer indicating the width of the window "
                "used for peak profiles. Peaks will be centered "
                "on their summits and include half of the window "
                "size upstream and half downstream of this point.",
                default=400,
            )
            shift_size = StringField(
                label="Shift size",
                description="Vector of values to try when computing optimal "
                "shift sizes. It should be specifeird as "
                "consecutive numbers vector with start:end",
                default="1:300",
            )

        advanced = GroupField(
            Advanced,
            label="Advanced parameters",
        )

    class Output:
        """Output fields to process ChipQC."""

        report_folder = DirField(label="ChipQC report folder")
        ccplot = FileField(label="Cross coverage score plot")
        coverage_histogram = FileField(label="SSD metric plot")
        peak_profile = FileField(label="Peak profile plot")
        peaks_barplot = FileField(label="Barplot of reads in peaks")
        peaks_density_plot = FileField(label="Density plot of reads in peaks")
        enrichment_heatmap = FileField(
            label="Heatmap of reads in genomic features", required=False
        )
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if inputs.alignment.entity_name != inputs.peaks.entity_name:
            self.error(
                "Sample names of aligned reads and called peaks do not "
                f" match. Alignment has {inputs.alignment.entity_name}, "
                f"while called peaks have {inputs.peaks.entity_name}."
            )

        basename = Path(inputs.alignment.bam.path).name
        assert basename.endswith(".bam")
        name = basename[:-4]
        report_folder = f"{name}_ChIPQCreport"

        if inputs.peaks.type == "data:chipseq:callpeak:macs14:":
            peaks_basename = Path(inputs.peaks.peaks_bed.path).name
            assert peaks_basename.endswith(".bed.gz")
            peaks_name = peaks_basename[:-7]
            peaks_file = f"{peaks_name}.bed"
            (Cmd["pigz"]["-cd", inputs.peaks.peaks_bed.path] > peaks_file)()
        elif inputs.peaks.type == "data:chipseq:callpeak:macs2:":
            peaks_file = inputs.peaks.called_peaks.path

        if inputs.alignment.build != inputs.peaks.build:
            self.error(
                "All input files must share the same build. "
                f"Aligment is based on {inputs.alignment.build} while "
                f"called peaks are based on {inputs.peaks.build}."
            )

        build = inputs.alignment.build
        genome_list = ["hg19", "hg38", "hg18", "mm10", "mm9", "rn4", "ce6", "dm3"]
        annotation = "NULL"
        if inputs.calculate_enrichment:
            if build in genome_list:
                annotation = f'"{build}"'
            else:
                self.warning(
                    f"Annotation for {build} is not supported. "
                    f'Supported builds are: {", ".join(genome_list)}.'
                )

        if inputs.blacklist:
            blacklist = f'"{inputs.blacklist.bed.path}"'
        else:
            blacklist = "NULL"

        r_input = (
            'library("ChIPQC"); '
            "sample <- ChIPQCsample("
            f'"{inputs.alignment.bam.path}",'
            f'peaks="{peaks_file}",'
            f"annotation={annotation},"
            f"mapQCth={inputs.advanced.quality_threshold},"
            f"blacklist={blacklist},"
            f"profileWin={inputs.advanced.profile_window},"
            f"shifts={inputs.advanced.shift_size});"
            "ChIPQCreport("
            "sample,"
            f'reportName="{name}_ChIPQC",'
            f'reportFolder="{report_folder}")'
        )

        return_code, _, _ = Cmd["Rscript"]["-e"][r_input] & TEE(retcode=None)

        if return_code:
            self.error("Error while running ChIPQC.")

        plots = Path(report_folder).glob("*.png")
        for plot_path in plots:
            copyfile(plot_path, f"{plot_path.stem}_mqc.png")

        outputs.ccplot = "CCPlot_mqc.png"
        outputs.coverage_histogram = "CoverageHistogramPlot_mqc.png"
        outputs.peak_profile = "PeakProfile_mqc.png"
        outputs.peaks_barplot = "Rip_mqc.png"
        outputs.peaks_density_plot = "Rap_mqc.png"
        outputs.report_folder = report_folder
        outputs.species = inputs.alignment.species
        outputs.build = build
        if Path("GenomicFeatureEnrichment_mqc.png").is_file():
            outputs.enrichment_heatmap = "GenomicFeatureEnrichment_mqc.png"

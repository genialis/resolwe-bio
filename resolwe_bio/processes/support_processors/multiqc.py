"""Prepare MultiQC report."""
import json
import os

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    DirField,
    FileHtmlField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    StringField,
)


def create_symlink(src, dst):
    """Create a symbolic link."""
    return Cmd["ln"]("-s", "--backup=numbered", src, dst)


def create_summary_table(samples, species, build):
    """Prepare sample summary MultiQC table."""
    sample_summary_json = {
        "id": "sample_info",
        "section_name": "Sample Info",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, sample_species, sample_build in zip(samples, species, build):
        if sample_build not in ["rRNA", "globin"]:
            sample_summary_json["data"][sample_name] = {
                "Species": sample_species,
                "Genome Build": sample_build,
            }

    with open("sample_data_mqc.json", "w") as out_file:
        json.dump(sample_summary_json, out_file)


def parse_chip_qc_report(report):
    """Parse ChiP-seq QC report file."""
    df = pd.read_csv(report, sep="\t")
    df.fillna("", inplace=True)
    return df.to_dict(orient="records")[0]


def create_prepeak_table(sample_names, reports):
    """Prepare ChIP-seq pre-peak MultiQC table."""
    prepeak_qc_json = {
        "pconfig": {"format": "{:,.2f}"},
        "id": "chip_seq_prepeak_qc",
        "section_name": "ChIP-seq pre-peak QC",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_chip_qc_report(report)
        prepeak_qc_json["data"][sample_name] = report_data

    with open("chipseq_prepeak_qc_mqc.json", "w") as out_file:
        json.dump(prepeak_qc_json, out_file)


def create_postpeak_table(sample_names, reports):
    """Prepare ChIP-seq pre-peak MultiQC table."""
    postpeak_qc_json = {
        "pconfig": {"format": "{:,.2f}"},
        "id": "chip_seq_postpeak_qc",
        "section_name": "ChIP-seq post-peak QC",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_chip_qc_report(report)
        postpeak_qc_json["data"][sample_name] = report_data

    with open("chipseq_postpeak_qc_mqc.json", "w") as out_file:
        json.dump(postpeak_qc_json, out_file)


def create_lib_strand_table(samples, reports):
    """Prepare library strandedness MultiQC table."""
    strand_codes = {
        "IU": "Strand non-specific (paired-end; -fr-unstranded)",
        "U": "Strand non-specific (single-end; -fr-unstranded)",
        "ISF": "Strand-specific forward (paired-end; -fr-secondstrand)",
        "OSF": "Strand-specific forward (paired-end; outward facing reads)",
        "SF": "Strand-specific forward (single-end; -fr-secondstrand)",
        "ISR": "Strand-specific reverse (paired-end; -fr-firststrand)",
        "OSR": "Strand-specific reverse (paired-end; outward facing reads)",
        "SR": "Strand-specific reverse (single-end; -fr-firststrand)",
    }

    lib_strand_json = {
        "id": "lib_strandedness",
        "section_name": "Library Strandedness",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(samples, reports):
        with open(report) as infile:
            data = json.load(infile)
            if "expected_format" in data:
                strandedness = data["expected_format"]
            else:
                raise ValueError("Cannot parse library type information file.")

        lib_strand_json["data"][sample_name] = {
            "Strandedness code": strandedness,
            "Description": strand_codes[strandedness],
        }

    with open("lib_strandedness_mqc.json", "w") as out_file:
        json.dump(lib_strand_json, out_file)


def process_strand_report_file(data, lib_type_samples, lib_type_reports):
    """Process Strandedness report file if it exists as Data output file."""
    try:
        if os.path.isfile(data.output.strandedness_report.path):
            lib_type_samples.append(data.entity.name)
            lib_type_reports.append(data.output.strandedness_report.path)
    except AttributeError:
        pass


def parse_bsrate_report(report):
    """Parse bsrate report file."""
    bsrate_dict = {}
    with open(report, "r") as r:
        first_line = r.readline()
        if not first_line.startswith("Bisulfite conversion rate process skipped."):
            bsrate_dict["Overall conversion rate"] = first_line.split()[-1]
            for strand in ["positive", "negative"]:
                line = r.readline().split()
                bsrate_dict[f"Conversion rate on {strand}"] = line[-2]
                bsrate_dict[f"Number of nucleotides used on {strand} strand"] = line[-1]
    return bsrate_dict


def create_bsrate_table(samples, reports):
    """Prepare bisulfite sequencing conversion rate MultiQC table."""
    bsrate_json = {
        "id": "wgbs_bsrate",
        "section_name": "WGBS conversion rate ",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }
    for sample_name, report in zip(samples, reports):
        report_data = parse_bsrate_report(report)
        if report_data:
            bsrate_json["data"][sample_name] = report_data
    if bsrate_json["data"]:
        with open("wgbs_bsrate_mqc.json", "w") as out_file:
            json.dump(bsrate_json, out_file)


def process_markdup_file(report):
    """Process samtools markdup file."""
    df = pd.read_csv(
        report, skiprows=0, sep=": ", header=0, engine="python"
    ).transpose()
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    data_table = df.to_dict(orient="records")[0]

    metrics = {}
    metrics["UNIQUE PAIRS"] = data_table["PAIRED"] - data_table["DUPLICATE PAIR"]
    metrics["UNIQUE SINGLE"] = data_table["SINGLE"] - data_table["DUPLICATE SINGLE"]
    metrics["DUPLICATE PAIRS OPTICAL"] = data_table["DUPLICATE PAIR OPTICAL"]
    metrics["DUPLICATE PAIRS NONOPTICAL"] = (
        data_table["DUPLICATE PAIR"] - data_table["DUPLICATE PAIR OPTICAL"]
    )
    metrics["DUPLICATE SINGLE OPTICAL"] = data_table["DUPLICATE SINGLE OPTICAL"]
    metrics["DUPLICATE SINGLE NONOPTICAL"] = (
        data_table["DUPLICATE SINGLE"] - data_table["DUPLICATE SINGLE OPTICAL"]
    )
    metrics["EXCLUDED"] = data_table["EXCLUDED"]
    return metrics


def create_markdup_plot(samples, reports):
    """Prepare samtools markdup MultiQC table."""
    markdup_json = {
        "id": "samtools_markdup",
        "section_name": "Samtools markdup statistics",
        "description": "*Please note that excluded reads are those marked as secondary, supplementary, QC failed or "
        "unmapped reads and are not used for calculating duplicates.",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "rmdup_bargraph",
            "title": "Samtools deduplication stats",
            "ylab": "Number of reads",
        },
        "file_format": "json",
        "data": {},
    }
    for sample_name, report in zip(samples, reports):
        report_data = process_markdup_file(report)
        if report_data:
            markdup_json["data"][sample_name] = report_data
    if markdup_json["data"]:
        with open("wgbs_markdup_mqc.json", "w") as out_file:
            json.dump(markdup_json, out_file)


def parse_nanostring_report(report):
    """Parse Nanostring sample QC report file."""
    df = pd.read_csv(
        report,
        sep="\t",
        names=[
            "Samples",
            "Zero Counts",
            "Signal",
            "Neg. Controls Mean",
            "Neg. Controls Std Dev",
            "Noise Cutoff",
        ],
        usecols=[
            "Zero Counts",
            "Signal",
            "Neg. Controls Mean",
            "Neg. Controls Std Dev",
            "Noise Cutoff",
        ],
        header=0,
    )
    df.fillna("NA", inplace=True)
    return df.to_dict(orient="records")[0]


def create_nanostring_table(sample_names, reports):
    """Prepare Nanostring Sample QC MultiQC table."""
    sample_qc_json = {
        "pconfig": {"format": "{:,.2f}"},
        "id": "nanostring_sample_qc",
        "section_name": "Nanostring sample QC",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_nanostring_report(report)
        sample_qc_json["data"][sample_name] = report_data

    with open("nanostring_sample_qc_mqc.json", "w") as out_file:
        json.dump(sample_qc_json, out_file)


def parse_lane_attributes(report):
    """Parse Nanostring lane attributes report file."""
    df = pd.read_csv(report, sep="\t", index_col=0, header=None).T
    df.drop(df.columns[0], axis=1, inplace=True)
    return df.to_dict(orient="records")[0]


def create_lane_table(sample_names, reports):
    """Prepare Nanostring lane attributes MultiQC table."""
    lane_json = {
        "pconfig": {"format": "{}", "scale": False},
        "id": "nanostring_lane_attributes",
        "section_name": "Nanostring lane attributes",
        "plot_type": "table",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(sample_names, reports):
        report_data = parse_lane_attributes(report)
        lane_json["data"][sample_name] = report_data

    with open("nanostring_lane_attributes_mqc.json", "w") as out_file:
        json.dump(lane_json, out_file)


class MultiQC(Process):
    """Aggregate results from bioinformatics analyses across many samples into a single report.

    [MultiQC](http://www.multiqc.info) searches a given directory for analysis logs and compiles a
    HTML report. It's a general purpose tool, perfect for summarising the output from numerous
    bioinformatics tools.
    """

    slug = "multiqc"
    process_type = "data:multiqc"
    name = "MultiQC"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.3.1"},
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "Other"
    data_name = "MultiQC report"
    version = "1.11.2"

    class Input:
        """Input fields to process MultiQC."""

        data = ListField(
            DataField(
                data_type="",
                description="Select multiple data objects for which the MultiQC report is to be "
                "generated.",
            ),
            label="Input data",
        )

        class Advanced:
            """Options."""

            dirs = BooleanField(
                label="--dirs",
                default=True,
                description="Prepend directory to sample names.",
            )

            dirs_depth = IntegerField(
                label="--dirs-depth",
                default=-1,
                description="Prepend a specified number of directories to sample names. Enter a "
                "negative number (default) to take from start of path.",
            )

            fullnames = BooleanField(
                label="--fullnames",
                default=False,
                description="Disable the sample name cleaning (leave as full file name).",
            )

            config = BooleanField(
                label="Use configuration file",
                default=True,
                description="Use Genialis configuration file for MultiQC report.",
            )

            cl_config = StringField(
                label="--cl-config",
                required=False,
                description="Enter text with command-line configuration options to override the "
                "defaults (e.g. custom_logo_url: https://www.genialis.com).",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields."""

        report = FileHtmlField(label="MultiQC report")
        report_data = DirField(label="Report data")

    def run(self, inputs, outputs):
        """Run the analysis."""
        samples = []
        species = []
        build = []
        lib_type_samples = []
        lib_type_reports = []
        chip_seq_samples = []
        chip_seq_prepeak_reports = []
        chip_seq_postpeak_samples = []
        chip_seq_postpeak_reports = []
        bsrate_samples = []
        bsrate_reports = []
        markdup_samples = []
        markdup_reports = []
        rcc_samples = []
        rcc_reports = []
        rcc_lane_reports = []
        unsupported_data = []

        for d in inputs.data:
            try:
                sample_dir = d.entity.name
                os.makedirs(sample_dir, exist_ok=True)
                if d.entity.name and d.output.species and d.output.build:
                    samples.append(d.entity.name)
                    species.append(d.output.species)
                    build.append(d.output.build)
            except AttributeError:
                pass

            if d.process.type.startswith("data:reads:fastq:single"):
                for fq_report in d.output.fastqc_archive:
                    name = os.path.basename(fq_report.path)
                    create_symlink(fq_report.path, os.path.join(sample_dir, name))

            elif d.process.type.startswith("data:reads:fastq:paired"):
                for fq_report in d.output.fastqc_archive + d.output.fastqc_archive2:
                    name = os.path.basename(fq_report.path)
                    create_symlink(fq_report.path, os.path.join(sample_dir, name))

            elif d.process.type.startswith("data:alignment:bam:markduplicate"):
                name = os.path.basename(d.output.metrics_file.path)
                create_symlink(
                    d.output.metrics_file.path, os.path.join(sample_dir, name)
                )

            elif d.process.type == "data:alignment:bam:star:":
                stats_file = os.path.basename(d.output.stats.path)
                assert stats_file.endswith("_stats.txt")
                bam_name = stats_file[:-10]

                if d.output.build == "rRNA":
                    rrna_report = f"{bam_name}.rRNA.Log.final.out"
                    create_symlink(
                        d.output.stats.path, os.path.join(sample_dir, rrna_report)
                    )
                elif d.output.build == "globin":
                    globin_report = f"{bam_name}.globin.Log.final.out"
                    create_symlink(
                        d.output.stats.path, os.path.join(sample_dir, globin_report)
                    )
                else:
                    report = f"{bam_name}.Log.final.out"
                    create_symlink(
                        d.output.stats.path, os.path.join(sample_dir, report)
                    )

            elif d.process.type == "data:alignment:bam:walt:":
                try:
                    if os.path.isfile(d.output.duplicates_report.path):
                        dup_report_path = d.output.duplicates_report.path
                        name = os.path.basename(dup_report_path)
                        create_symlink(dup_report_path, os.path.join(sample_dir, name))
                        markdup_samples.append(d.entity.name)
                        markdup_reports.append(dup_report_path)
                        create_markdup_plot(markdup_samples, markdup_reports)
                except AttributeError:
                    pass

            elif d.process.type.startswith("data:alignment:bam"):
                name = os.path.basename(d.output.stats.path)
                create_symlink(d.output.stats.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:expression:featurecounts:":
                name = os.path.basename(d.output.counts_summary.path)
                create_symlink(
                    d.output.counts_summary.path, os.path.join(sample_dir, name)
                )
                # Strandedness report exists only if auto detection was enabled
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.process.type == "data:chipseq:callpeak:macs2:":
                name = os.path.basename(d.output.called_peaks.path)
                create_symlink(
                    d.output.called_peaks.path, os.path.join(sample_dir, name)
                )
                chip_seq_samples.append(d.entity.name)
                chip_seq_prepeak_reports.append(d.output.case_prepeak_qc.path)
                chip_seq_postpeak_samples.append(d.entity.name)
                chip_seq_postpeak_reports.append(d.output.chip_qc.path)
                # MACS2 analysis can be run without the background sample,
                # thus the associated report might not exits
                try:
                    if os.path.isfile(d.output.control_prepeak_qc.path):
                        chip_seq_samples.append(f"Background of {d.entity.name}")
                        chip_seq_prepeak_reports.append(
                            d.output.control_prepeak_qc.path
                        )
                except AttributeError:
                    pass

            elif d.process.type == "data:samtools:idxstats:":
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:qorts:qc:":
                qorts_path = os.path.join(sample_dir, "QoRTs")
                os.mkdir(qorts_path)
                name = os.path.basename(d.output.summary.path)
                create_symlink(d.output.summary.path, os.path.join(qorts_path, name))

            elif d.process.type == "data:expression:salmon:":
                create_symlink(d.output.salmon_output.path, sample_dir)
                # Strandedness report might not exist in legacy Salmon objects
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.process.type.startswith("data:alleyoop"):
                create_symlink(d.output.report.path, sample_dir)
                # Alleyoop summary may contain PCA plot data
                try:
                    if os.path.isfile(d.output.plot_data.path):
                        create_symlink(d.output.plot_data.path, sample_dir)
                except AttributeError:
                    pass

            elif d.process.type.startswith("data:picard"):
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:wgbs:bsrate:":
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))
                bsrate_samples.append(d.entity.name)
                bsrate_reports.append(d.output.report.path)
                create_bsrate_table(bsrate_samples, bsrate_reports)

            elif d.process.type == "data:chipqc:":
                plot_paths = [
                    d.output.ccplot.path,
                    d.output.coverage_histogram.path,
                    d.output.peak_profile.path,
                    d.output.peaks_barplot.path,
                    d.output.peaks_density_plot.path,
                ]
                for path in plot_paths:
                    name = os.path.basename(path)
                    create_symlink(path, os.path.join(sample_dir, name))
                # ChipQC may contain enrichment heatmap
                try:
                    if os.path.isfile(d.output.enrichment_heatmap.path):
                        name = os.path.basename(d.output.enrichment_heatmap.path)
                        create_symlink(
                            d.output.enrichment_heatmap.path,
                            os.path.join(sample_dir, name),
                        )
                except AttributeError:
                    pass

            elif d.process.type == "data:nanostring:rcc:":
                # Sample_qc is an optional field
                try:
                    name = os.path.basename(d.output.sample_qc.path)
                    create_symlink(
                        d.output.sample_qc.path, os.path.join(sample_dir, name)
                    )
                    rcc_samples.append(d.entity.name)
                    rcc_reports.append(d.output.sample_qc.path)
                    create_nanostring_table(rcc_samples, rcc_reports)

                    lane_name = os.path.basename(d.output.lane_attributes.path)
                    create_symlink(
                        d.output.lane_attributes.path,
                        os.path.join(sample_dir, lane_name),
                    )
                    rcc_lane_reports.append(d.output.lane_attributes.path)
                    create_lane_table(rcc_samples, rcc_lane_reports)

                except AttributeError:
                    pass
            else:
                unsupported_data.append(d.name)

        if unsupported_data:
            ext = ", ..." if len(unsupported_data) > 5 else ""
            self.warning(
                f"The Input data {', '.join(unsupported_data[:5])}{ext} is not supported "
                f"by the MultiQC analysis."
            )

        create_summary_table(samples, species, build)

        if lib_type_samples and lib_type_reports:
            create_lib_strand_table(lib_type_samples, lib_type_reports)

        if chip_seq_samples and chip_seq_prepeak_reports:
            create_prepeak_table(chip_seq_samples, chip_seq_prepeak_reports)

        if chip_seq_postpeak_samples and chip_seq_postpeak_reports:
            create_postpeak_table(chip_seq_postpeak_samples, chip_seq_postpeak_reports)

        args = [
            "-dd",
            inputs.advanced.dirs_depth,
        ]

        if inputs.advanced.dirs:
            args.append("-d")

        if inputs.advanced.fullnames:
            args.append("-s")

        if inputs.advanced.config:
            args.extend(["-c", "/opt/resolwebio/assets/multiqc_config.yml"])

        if inputs.advanced.cl_config:
            args.extend(["--cl-config ", inputs.advanced.output.cl_config])

        with Cmd.env(LC_ALL="C.UTF-8"):
            return_code, _, _ = Cmd["multiqc"]["."][args] & TEE(retcode=None)
            if return_code:
                self.error("MultiQC analysis failed.")

        if not os.path.isdir("multiqc_data") and not os.path.isfile(
            "multiqc_report.html"
        ):
            self.error("MultiQC finished without creating outputs.")

        outputs.report = "multiqc_report.html"
        outputs.report_data = "multiqc_data"

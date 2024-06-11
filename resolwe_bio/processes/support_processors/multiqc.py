"""Prepare MultiQC report."""

import gzip
import json
import os
import shutil
import zipfile
from pathlib import Path

import pandas as pd
import yaml
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


def clean_name(sample_name, to_remove, error):
    """Clean a sample name."""
    for substr in to_remove:
        sample_name = sample_name.replace(substr, "")

    if not sample_name:
        error(
            "Sample name only contains elements which are removed during sample name cleanup."
            f"Avoid naming samples with just {' ,'.join(to_remove)}."
        )

    return sample_name


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


def parse_genebody_report(report):
    """Parse QoRTs gene body coverage metrics report file."""
    df = pd.read_csv(report, sep="\t", compression="gzip")
    df["QUANTILE"] *= 100
    dict = {k: v for k, v in zip(df["QUANTILE"], df["TOTAL"])}
    return dict


def create_coverage_plot(sample_names, reports):
    """Prepare QoRTs gene body coverage plot."""
    genebody_qc_json = {
        "id": "genebody_qc",
        "section_name": "QoRTs QC - gene body coverage information",
        "plot_type": "linegraph",
        "file_format": "json",
        "pconfig": {
            "ylab": "Proportion of Read-Pairs",
            "xlab": "Percentile of Gene Body (5' -> 3')",
        },
        "data": {},
    }
    source_fn = "qorts_output/QC.geneBodyCoverage.byExpr.avgPct.txt.gz"

    for sample_name, report in zip(sample_names, reports):
        with zipfile.ZipFile(report, "r") as archive:
            archive.extract(source_fn)

        report_data = parse_genebody_report(source_fn)
        genebody_qc_json["data"][sample_name] = report_data

    with open("qorts_genebody_mqc.json", "w") as out_file:
        json.dump(genebody_qc_json, out_file)


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


def sum_featurecounts_columns(summary_file, out_file):
    """Prepare input for featureCounts."""
    exp = pd.read_csv(
        summary_file,
        sep="\t",
        index_col="Status",
        dtype={
            "Status": str,
        },
    )
    if len(exp.columns) > 1:
        sum_column = exp.columns[0].split(":")[0]
        exp[sum_column] = exp.sum(axis=1)

        exp = exp.astype({sum_column: int})
        exp[[sum_column]].to_csv(
            out_file,
            index_label="Status",
            sep="\t",
        )


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


def parse_counts_summary(report):
    """Parse STAR quantification summary report file."""
    df = pd.read_csv(
        report,
        sep="\t",
        index_col=0,
    )
    assigned_reads = int(df.loc[df.index == "N_assigned", "Read count"].values)
    mapped_reads_sum = df.loc[
        ["N_multimapping", "N_noFeature", "N_ambiguous", "N_assigned"], "Read count"
    ].sum()
    mapped_reads_sum = float(mapped_reads_sum)
    percent_assigned = assigned_reads / mapped_reads_sum * 100

    out_dict = {
        "Assigned reads": assigned_reads,
        "% of assigned reads": percent_assigned,
    }

    return out_dict


def update_generalstats_table(sample_names, reports):
    """Update general statistics table with new information."""
    counts_json = {
        "section_name": "STAR quantification",
        "plot_type": "generalstats",
        "file_format": "json",
        "data": {},
    }

    for sample_name, report in zip(sample_names, reports):
        sample_name += f" | {sample_name}"
        report_data = parse_counts_summary(report)
        counts_stats = list(report_data)

        counts_json["data"][sample_name] = {
            k: report_data[k] for k in counts_stats if k in report_data
        }

    with open("STAR quantification_mqc.json", "w") as out_file:
        json.dump(counts_json, out_file)


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
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"},
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "QC"
    data_name = "MultiQC report"
    version = "1.25.0"

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
                default=-5,
                description="Prepend a specified number of directories to sample names. Enter a "
                "negative number to take from start of path.",
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
        star_quantification_samples = []
        star_quantification_reports = []
        qorts_samples = []
        qorts_reports = []

        config_file = "/opt/resolwebio/assets/multiqc_config.yml"
        with open(config_file) as handle:
            mqc_config = yaml.safe_load(handle)

        for d in inputs.data:
            try:
                # Here we have to remove some filename suffixes
                # to avoid missing data in the final report. This
                # workaround is used to avoid duplicated names caused by
                # file name cleaning.
                # For example, `some_sample.fastq.gz/stats_L001.txt`
                # and `some_sample.fastq.gz/stats_L002.txt` would
                # both be simplified to `some_sample` and only the first
                # would be included in the report.
                sample_name = clean_name(
                    sample_name=d.entity.name,
                    to_remove=mqc_config["extra_fn_clean_exts"],
                    error=self.error,
                )
                sample_dir = sample_name
                os.makedirs(sample_dir, exist_ok=True)
                if sample_name and d.output.species and d.output.build:
                    samples.append(sample_name)
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
                    count_report = "ReadsPerGene.out.tab.gz"
                    report = f"{bam_name}.Log.final.out"
                    alignment_count = 1
                    alignment_dir = os.path.join(sample_dir, f"STAR_{alignment_count}")
                    if not os.path.isdir(alignment_dir):
                        os.makedirs(alignment_dir)
                    else:
                        while os.path.isdir(alignment_dir):
                            alignment_count += 1
                            alignment_dir = os.path.join(
                                sample_dir, f"STAR_{alignment_count}"
                            )
                    os.makedirs(alignment_dir, exist_ok=True)
                    dst = os.path.join(alignment_dir, report)
                    create_symlink(d.output.stats.path, dst)

                if d.output.gene_counts:
                    outfile = os.path.join(alignment_dir, count_report)
                    if not os.path.isfile(outfile):
                        with gzip.open(d.output.gene_counts.path, "rb") as f_in:
                            with open(
                                os.path.join(alignment_dir, "ReadsPerGene.out.tab"),
                                "wb",
                            ) as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    else:
                        self.warning(f"File {outfile} already exists. Skipping.")

            elif d.process.type == "data:alignment:bam:walt:":
                try:
                    if os.path.isfile(d.output.duplicates_report.path):
                        dup_report_path = d.output.duplicates_report.path
                        name = os.path.basename(dup_report_path)
                        create_symlink(dup_report_path, os.path.join(sample_dir, name))
                        markdup_samples.append(sample_name)
                        markdup_reports.append(dup_report_path)
                        create_markdup_plot(markdup_samples, markdup_reports)
                except AttributeError:
                    pass

            elif d.process.type == "data:alignment:bam:bqsr:":
                name = os.path.basename(d.output.recal_table.path)
                create_symlink(
                    d.output.recal_table.path, os.path.join(sample_dir, name)
                )

            elif d.process.type.startswith("data:alignment:bam"):
                name = os.path.basename(d.output.stats.path)
                create_symlink(d.output.stats.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:expression:featurecounts:":
                name = os.path.basename(d.output.counts_summary.path)
                sum_featurecounts_columns(
                    summary_file=d.output.counts_summary.path,
                    out_file=os.path.join(sample_dir, f"summed_{name}"),
                )
                if not os.path.exists(os.path.join(sample_dir, f"summed_{name}")):
                    create_symlink(
                        d.output.counts_summary.path, os.path.join(sample_dir, name)
                    )
                # Strandedness report exists only if auto detection was enabled
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.process.type == "data:expression:star:":
                name = os.path.basename(d.output.counts_summary.path)
                create_symlink(d.output.counts_summary, os.path.join(sample_dir, name))
                star_quantification_samples.append(sample_name)
                star_quantification_reports.append(d.output.counts_summary.path)
                update_generalstats_table(
                    sample_names=star_quantification_samples,
                    reports=star_quantification_reports,
                )

            elif d.process.type == "data:chipseq:callpeak:macs2:":
                name = os.path.basename(d.output.called_peaks.path)
                create_symlink(
                    d.output.called_peaks.path, os.path.join(sample_dir, name)
                )
                chip_seq_samples.append(sample_name)
                chip_seq_prepeak_reports.append(d.output.case_prepeak_qc.path)
                # When MACS2 analysis is run in broad peak mode (--broad), the postpeak
                # reports are not generated
                try:
                    if os.path.isfile(d.output.chip_qc.path):
                        chip_seq_postpeak_samples.append(sample_name)
                        chip_seq_postpeak_reports.append(d.output.chip_qc.path)
                except AttributeError:
                    pass
                # MACS2 analysis can be run without the background sample,
                # thus the associated report might not exits
                try:
                    if os.path.isfile(d.output.control_prepeak_qc.path):
                        chip_seq_samples.append(f"Background of {sample_name}")
                        chip_seq_prepeak_reports.append(
                            d.output.control_prepeak_qc.path
                        )
                except AttributeError:
                    pass

            elif d.process.type == "data:samtools:idxstats:":
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:qorts:qc:":
                summary_name = os.path.basename(d.output.summary.path)
                create_symlink(
                    d.output.summary.path, os.path.join(sample_dir, summary_name)
                )

                report_path = d.output.qorts_data.path
                qorts_samples.append(sample_name)
                qorts_reports.append(report_path)
                create_symlink(
                    report_path, os.path.join(sample_dir, os.path.basename(report_path))
                )

            elif d.process.type == "data:rnaseqc:qc:":
                name = os.path.basename(d.output.metrics.path)
                create_symlink(
                    src=d.output.metrics.path, dst=os.path.join(sample_dir, name)
                )

            elif d.process.type == "data:expression:salmon:":
                # Symlink files/dirs without the parent directory to
                # attach it to the same sample in the general summary.
                for out_file in Path(d.output.salmon_output.path).iterdir():
                    create_symlink(str(out_file), str(Path(sample_dir) / out_file.name))
                # Strandedness report might not exist in legacy Salmon objects
                process_strand_report_file(d, lib_type_samples, lib_type_reports)

            elif d.process.type.startswith("data:picard"):
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))

            elif d.process.type == "data:wgbs:bsrate:":
                name = os.path.basename(d.output.report.path)
                create_symlink(d.output.report.path, os.path.join(sample_dir, name))
                bsrate_samples.append(sample_name)
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
                            src=d.output.enrichment_heatmap.path,
                            dst=os.path.join(sample_dir, name),
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

        if qorts_samples and qorts_reports:
            create_coverage_plot(sample_names=qorts_samples, reports=qorts_reports)

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
            args.extend(["-c", config_file])

        if inputs.advanced.cl_config:
            args.extend(["--cl-config", inputs.advanced.cl_config])

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

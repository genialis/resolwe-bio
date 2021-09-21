"""Prepare QC report."""
import os

import pandas as pd
import plotly.graph_objects as go

HTML_TEMPLATE = """
<html>
<head>
<title>mRNA QC report</title>
<meta charset="utf-8" />
</head>
<body>
<h1>mRNA QC report</h1>
{{CONTENT}}
</body>
</html>
"""

from resolwe.process import (
    DataField,
    DirField,
    FileField,
    FileHtmlField,
    GroupField,
    ListField,
    Process,
)


def summarize_reads_number(data, name):
    """Summarize reads for mates of a sample."""
    sum_reads = (
        data.loc[:, [c for c in data.columns if c.endswith("-total_sequences")]]
        .sum(axis=0)
        .rename(name)
        .T
    )
    return sum_reads


def average_mates(data, name):
    """Compute avregae values for mates of a sample."""
    mean = (
        data.drop(columns=[c for c in data.columns if c.endswith("-total_sequences")])
        .mean(axis=0)
        .rename(name)
        .T
    )
    return mean


def create_table(stats_sum, stats_mean):
    """Creating pandas table."""
    general_mean = pd.concat(stats_mean, axis=1).T
    general_mean.index.name = "Sample"
    general_sum = pd.concat(stats_sum, axis=1).T
    general_sum.index.name = "Sample"
    general_stats = pd.concat([general_sum, general_mean], axis=1)
    return general_stats


def filter_table(table):
    """Filtering pandas table."""
    fastqc = table.loc[:, table.columns.str.contains("FastQC")]
    star = table.loc[:, table.columns.str.contains("STAR")]
    if list(table.loc[:, table.columns.str.contains("Salmon")].columns) > []:
        quant = table.loc[:, table.columns.str.contains("Salmon")]
    elif list(table.loc[:, table.columns.str.contains("featureCounts")].columns) > []:
        quant = table.loc[:, table.columns.str.contains("featureCounts")]
    return fastqc, star, quant


def get_violin_figure(data, title, ytitle):
    """Get violin plot."""
    fig = go.Figure()
    for y, name in data:
        fig.add_trace(
            go.Violin(
                y=y,
                x0=name,
                name=name,
                hovertext=y.index,
                line_color="Black",
                fillcolor="#E5D2B4",
                box_fillcolor="gray",
                meanline_visible=True,
                box_visible=True,
                points="all",
                jitter=0.1,
            )
        )
    fig.update_layout(
        autosize=True,
        template="simple_white",
        title=title,
    )
    fig.update_yaxes(
        title=ytitle,
        showgrid=True,
        gridcolor="lightgrey",
    )

    return fig


def get_barplot_figure(data, title, ytitle, thr):
    """Get bar plot."""
    fig = go.Figure()
    for y, name, color in data:
        fig.add_trace(
            go.Bar(
                x=y.index,
                y=y,
                name=name,
                marker_color=color,
            )
        )
        fig.update_layout(
            autosize=True,
            template="simple_white",
            barmode="group",
            title=title,
            xaxis_title="Sample",
            yaxis_title=ytitle,
        )
    return fig


def create_report(path, HTML_TEMPLATE, html_content):
    """Write all to .html file"""
    with open(f"{path}/qc_report.html", "wt") as f:
        f.write(HTML_TEMPLATE.replace("{{CONTENT}}", html_content)).close()
    return f


class QCReport(Process):
    """Aggregate QC results from [MultiQC](http://www.multiqc.info) reports across many samples into a
    single HTML QC report.
    """

    slug = "qc_report"
    process_type = "data:qc"
    name = "QC Report"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:2.8.0"},
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
    data_name = "QC report"
    version = "1.0.0"
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields to generate a QC Report."""

        data = ListField(
            DataField(
                data_type="data:multiqc",
                description="Select MultiQC data objects for which the QC report is to be "
                "generated.",
            ),
            label="Input data",
        )

    class Output:
        """Output fields."""

        report = FileHtmlField(label="QC report")
        table = FileField(label="General statistics table")
        report_data = DirField(label="Report data")

    def run(self, inputs, outputs):
        """Run the analysis."""

        stats_mean = []
        stats_sum = []
        unsupported_data = []

        for d in inputs.data:
            # try:
            if d.process.type.startswith("data:multiqc"):
                # for d.output.multiqc_data:
                path = os.path.basename(d.output.multiqc_data.path)
                data = pd.read_csv(
                    f"{path}/multiqc_general_stats.txt",
                    sep="\t",
                    index_col="Sample",
                )
                name = d.name
                stats_sum.append(summarize_reads_number(data, name))
                stats_mean.append(average_mates(data, name))
            # except AttributeError:
            #     pass

            else:
                unsupported_data.append(d.name)

            if unsupported_data:
                ext = ", ..." if len(unsupported_data) > 5 else ""
            self.warning(
                f"The Input data {', '.join(unsupported_data[:5])}{ext} is not supported "
                f"by the QC report process."
            )

        # Creating general statistics table
        general_statistics = create_table(stats_sum, stats_mean)

        # Exporting to .tsv
        os.makedir(report_data)
        os.chdir(report_data)
        general_statistics.to_csv("general_stats_table.tsv", sep="\t")

        # Filter table
        fastqc, star, quant = filter_table(general_statistics)

        html_content = ""
        to_html_kwargs = {
            "full_html": False,
            "include_plotlyjs": "cdn",
        }

        # FASTQ - Number of raw and trimmed read counts
        data = (
            [
                (
                    fastqc.loc[:, fastqc.columns.str.contains("raw-total_sequences")],
                    "Raw reads",
                    "#E5D2B4",
                ),
                (
                    fastqc.loc[
                        :, fastqc.columns.str.contains("trimmed-total_sequences")
                    ],
                    "Trimmed reads",
                    "#A98A53",
                ),
            ],
        )
        fig = get_violin_figure(
            data, title="Total number of reads per sample", ytitle="Number of reads"
        )
        html_content += fig.to_html(**to_html_kwargs)
        if len(fastqc.index) <= 100:
            fig = get_bar_figure(
                data,
                title="Total number of reads per sample",
                ytitle="Number of reads",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # FASTQ - GC content
        data = (
            [
                (
                    fastqc.loc[:, fastqc.columns.str.contains("raw-percent_gc")],
                    "Raw reads",
                    "#E5D2B4",
                ),
                (
                    fastqc.loc[:, fastqc.columns.str.contains("trimmed-percent_gc")],
                    "Trimmed reads",
                    "#A98A53",
                ),
            ],
        )
        fig = get_violin_figure(
            data,
            title="GC content [%]",
            ytitle="GC content [%]",
        )
        html_content += fig.to_html(**to_html_kwargs)
        if len(fastqc.index) <= 100:
            fig = get_bar_figure(
                data,
                title="GC content [%]",
                ytitle="GC content [%]",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # FASTQ - Duplication
        data = (
            [
                (
                    fastqc.loc[
                        :, fastqc.columns.str.contains("raw-percent_duplicates")
                    ],
                    "Raw reads",
                    "#E5D2B4",
                ),
                (
                    fastqc.loc[
                        :, fastqc.columns.str.contains("trimmed-percent_duplicates")
                    ],
                    "Trimmed reads",
                    "#A98A53",
                ),
            ],
        )
        fig = get_violin_figure(
            data,
            title="Sequence duplication [%]",
            ytitle="Sequence duplication [%]",
        )
        html_content += fig.to_html(**to_html_kwargs)
        if len(fastqc.index) <= 100:
            fig = get_bar_figure(
                data,
                title="Sequence Duplication [%]",
                ytitle="Sequence duplication [%]",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # FASTQ - Sequence length
        data = (
            [
                (
                    fastqc.loc[
                        :, fastqc.columns.str.contains("raw-avg_sequence_length")
                    ],
                    "Raw reads",
                    "#E5D2B4",
                ),
                (
                    fastqc.loc[
                        :, fastqc.columns.str.contains("trimmed-avg_sequence_length")
                    ],
                    "Trimmed reads",
                    "#A98A53",
                ),
            ],
        )
        fig = get_violin_figure(
            data,
            title="Sequence length",
            ytitle="Sequence length",
        )
        html_content += fig.to_html(**to_html_kwargs)
        if len(fastqc.index) <= 100:
            fig = get_bar_figure(
                data,
                title="Sequence length",
                ytitle="Sequence length",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # Alignment - number uniquely mapped
        if quant.loc[:, quant.columns.str.contains("Salmon")]:
            data = (
                [
                    (
                        quant.loc[:, quant.columns.str.contains("salmon-num_mapped")],
                        "Uniquely mapped",
                        "#E5D2B4",
                    )
                ],
            )
            data_per = (
                [
                    (
                        quant.loc[
                            :, quant.columns.str.contains("salmon-percent_mapped")
                        ],
                        "Uniquely mapped",
                        "#E5D2B4",
                    )
                ],
            )

        elif quant.loc[:, quant.columns.str.contains("featureCounts")]:
            data = [
                (
                    star.loc[:, star.columns.str.contains("star-uniquely_mapped")],
                    "Uniquely mapped",
                    "#E5D2B4",
                )
            ]
            data_per = [
                (
                    star.loc[
                        :, star.columns.str.contains("star-uniquely_mapped_percent")
                    ],
                    "Uniquely mapped",
                    "#E5D2B4",
                )
            ]

        fig = get_violin_figure(
            data,
            title="Alignment statistics - number of mapped reads",
            ytitle="Number of mapped reads",
        )
        html_content += fig.to_html(**to_html_kwargs)

        fig = get_violin_figure(
            data_per,
            title="Alignment statistics - percent of mapped reads",
            ytitle="Percent of mapped reads",
        )
        html_content += fig.to_html(**to_html_kwargs)

        if len(quant.index) <= 100:
            fig = get_bar_figure(
                data,
                title="Alignment statistics - number of mapped reads",
                ytitle="Number of mapped reads",
            )
            html_content += fig.to_html(**to_html_kwargs)
            fig = get_bar_figure(
                data_per,
                title="Alignment statistics - percent of mapped reads",
                ytitle="Percent of mapped read",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # Alignment - percent uniquely mapped to rRNA, Globin:
        data = [
            (
                star.loc[:, star.columns.str.contains("rrna-uniquely_mapped_percent")],
                "rRNA",
                "#E5D2B4",
            ),
            (
                star.loc[
                    :, star.columns.str.contains("globin-uniquely_mapped_percent")
                ],
                "Globin",
                "#A98A53",
            ),
        ]
        fig = get_violin_figure(
            data,
            title="Alignment statistics - % of mapped reads to rRNA, Globin",
            ytitle="Percent of mapped reads",
        )
        html_content += fig.to_html(**to_html_kwargs)

        if len(star.index) <= 100:
            fig = get_bar_figure(
                data,
                title="Alignment statistics - % of mapped reads to rRNA, Globin",
                ytitle="Percent of mapped reads",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # Quantification
        if quant.loc[:, quant.columns.str.contains("featureCounts")]:
            data = (
                [
                    (
                        quant.loc[
                            :, quant.columns.str.contains("featurecounts-Assigned")
                        ],
                        "Assigned reads",
                        "#E5D2B4",
                    ),
                ],
            )
            data_per = (
                [
                    (
                        quant.loc[
                            :,
                            quant.columns.str.contains(
                                "featurecounts-percent_assigned"
                            ),
                        ],
                        "Assigned reads",
                        "#E5D2B4",
                    ),
                ],
            )
            fig = get_violin_figure(
                data,
                title="Quantification statistics - number of assigned reads",
                ytitle="Number of assigned reads",
            )
            html_content += fig.to_html(**to_html_kwargs)

            fig = get_violin_figure(
                data_per,
                title="Quantification statistics - % of assigned reads",
                ytitle="Percent of assigned reads",
            )
            html_content += fig.to_html(**to_html_kwargs)

            if len(quant.index) <= 100:
                fig = get_bar_figure(
                    data,
                    title="Quantification statistics - number of assigned reads",
                    ytitle="Number of assigned reads",
                )
                html_content += fig.to_html(**to_html_kwargs)
                fig = get_bar_figure(
                    data_per,
                    title="Quantification statistics - % of assigned reads",
                    ytitle="Percent of assigned reads",
                )
                html_content += fig.to_html(**to_html_kwargs)

        # Salmon - genes with nonzero counts & number of chromosomes covered
        if quant.loc[:, quant.columns.str.contains("Salmon")]:
            data = (
                [
                    (
                        quant.loc[
                            :, quant.columns.str.contains("PercentWithNonzeroCounts")
                        ],
                        "Percent of genes with nonzero counts",
                        "#E5D2B4",
                    ),
                ],
            )
            data_c = (
                [
                    (
                        quant.loc[
                            :, quant.columns.str.contains("NumberOfChromosomesCovereds")
                        ],
                        "Number of chromosomes covered",
                        "#E5D2B4",
                    ),
                ],
            )
            fig = get_violin_figure(
                title="Quantification statistics - percent of genes with non-zero counts",
                ytitle="Genes with non-zero count",
            )
            html_content += fig.to_html(**to_html_kwargs)
            fig = get_violin_figure(
                data_c,
                title="Quantification statistics - number of chromosomes covered",
                ytitle="Number of chromosomes covered",
            )
            html_content += fig.to_html(**to_html_kwargs)

        # Create QC report -> write to .html
        report = create_report(path, HTML_TEMPLATE, html_content)

        if (
            not os.path.isdir("qc_report_data")
            and not os.path.isfile("qc_report.html")
            and not os.path.isfile("general_statistics.tsv")
        ):
            self.error("QC report generation finished without creating outputs.")

        outputs.table = "general_stats_table.tsv"
        outputs.report = "qc_report.html"
        outputs.report_data = "qc_report_data"

#!/usr/bin/env python3
"""Merge ChIP/ATAC-seq prepeak and postpeak QC reports."""

import argparse
from collections import defaultdict

import pandas as pd

parser = argparse.ArgumentParser(description="Merge ChIP/ATAC-Seq QC reports.")

parser.add_argument(
    "-f",
    "--file_path",
    required=True,
    nargs="+",
    help="List with paths to QC report files.",
)
parser.add_argument(
    "-n", "--sample_names", required=True, nargs="+", help="List of sample names."
)
parser.add_argument(
    "-r", "--report_type", required=True, nargs="+", help="List of report types."
)


if __name__ == "__main__":

    args = parser.parse_args()

    data = defaultdict(dict)
    for report_file, report_type, sample_name in zip(
        args.file_path, args.report_type, args.sample_names
    ):
        data[sample_name][report_type] = report_file

    prepeak_list = []
    postpeak_list = []
    for sample in data:
        if "prepeak" in data[sample]:
            prepeak = pd.read_csv(data[sample]["prepeak"], sep="\t")
            prepeak.index = [sample]
            prepeak_list.append(prepeak)
        if "postpeak" in data[sample]:
            postpeak = pd.read_csv(data[sample]["postpeak"], sep="\t")
            postpeak.index = [sample]
            postpeak_list.append(postpeak)

        if prepeak_list and postpeak_list:
            prepeaks = pd.concat(prepeak_list)
            postpeaks = pd.concat(postpeak_list)
            report = pd.merge(
                prepeaks, postpeaks, left_index=True, right_index=True, how="outer"
            )
        elif prepeak_list:
            report = pd.concat(prepeak_list)
        else:
            report = pd.concat(postpeak_list)

    report.to_csv(
        "QC_report.txt",
        sep="\t",
        na_rep="N/A",
        index_label="SAMPLE_NAME",
        float_format="%.3f",
    )

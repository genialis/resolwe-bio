#!/usr/bin/env python3
"""
Correct bed file to be compatible with bedToBigBed tool.

    1.) restrict all scores to the maximal value of 1000,
    2.) in strand column replace '?' with '.'.
"""

import argparse

import pandas as pd
from pandas.errors import EmptyDataError
from resolwe_runtime_utils import error, send_message

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-f", "--bed_file", help="Bed file.")
args = parser.parse_args()

try:
    df = pd.read_csv(args.bed_file, delimiter="\t", header=None, dtype=str)
except EmptyDataError:
    send_message(
        error(
            f"The input BED file {args.bed_file} is empty. Your analysis might "
            f"have failed to identify regions of interest (peaks, junctions, etc.)."
        )
    )
else:
    df.iloc[:, 4] = pd.to_numeric(df.iloc[:, 4]).round().astype(int)
    df.iloc[:, 4] = df.iloc[:, 4].clip(upper=1000)

    # if strand column exist replace '?' with '.'
    if len(df.columns) >= 6:
        df.iloc[:, 5] = df.iloc[:, 5].replace("?", ".")

    output_name = "_".join(["corrected", args.bed_file])
    df.to_csv(output_name, sep="\t", index=False, header=False)

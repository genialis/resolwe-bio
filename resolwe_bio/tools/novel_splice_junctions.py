#!/usr/bin/env python3
"""Filter novel splice junctions in BED12 format."""

import argparse
import os
import sys

import pandas as pd
from resolwe_runtime_utils import send_message, warning

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
    "-bed", "--bed_file", required=True, help="All splice junctions in BED12 format"
)
parser.add_argument(
    "-sj", "--novel_sj", required=True, help="Table of annotated novel splice junctions"
)

if __name__ == "__main__":
    args = parser.parse_args()
    bed_file = args.bed_file

    if os.path.getsize(bed_file) == 0:
        send_message(warning("Bed file has no entries."))
        os.rename(bed_file, "novel_sj.bed")
        sys.exit(0)

    bed = pd.read_csv(args.bed_file, delimiter="\t", header=None, dtype=str)
    novel_sj = pd.read_csv(args.novel_sj, delimiter="\t", dtype=str)
    bed_novel_sj = bed[bed[3].isin(novel_sj["name"])]

    bed_novel_sj.to_csv("novel_sj.bed", sep="\t", index=False, header=False)

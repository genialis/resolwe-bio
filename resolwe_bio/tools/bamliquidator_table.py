#!/usr/bin/env python2
"""Rename column names in bamliquidator matrix table to show sample names istead of bam names."""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse

import pandas as pd

parser = argparse.ArgumentParser(
    description="Rename column names in bamliquidator matrix table."
)

parser.add_argument(
    "-b", "--bam_files", nargs="+", required=True, help="List of bam files"
)
parser.add_argument(
    "-s", "--sample_names", nargs="+", required=True, help="List of sample names"
)
parser.add_argument("-m", "--matrix_table", required=True, help="Matrix table")

args = parser.parse_args()

df = pd.read_table(args.matrix_table, sep="\t", index_col="GENE_ID")

keys = [("bin_1_" + x) for x in args.bam_files]
values = args.sample_names
dictionary = dict(zip(keys, values))

df.rename(columns=dictionary, inplace=True)

df.to_csv("table_output.csv", index=True, header=True, sep=str("\t"))

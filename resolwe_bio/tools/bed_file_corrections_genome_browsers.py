#!/usr/bin/env python3
"""Restrict all scores to the maximal value of 1000."""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Restrict all scores to the maximal value of 1000.")
parser.add_argument('-f', '--bed_file', help="Bed file.")
args = parser.parse_args()

df = pd.read_csv(args.bed_file, delimiter='\t', header=None, dtype=str)
df.iloc[:, 4] = pd.to_numeric(df.iloc[:, 4]).round().astype(int)
df.iloc[:, 4] = df.iloc[:, 4].clip(upper=1000)
output_name = '_'.join(['corrected', args.bed_file])
df.to_csv(output_name, sep='\t', index=False, header=False)

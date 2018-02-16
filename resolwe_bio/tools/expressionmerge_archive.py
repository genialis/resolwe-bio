#!/usr/bin/env python3
"""Merge columns of multiple experiments by gene id."""

import argparse
import csv
import os

import pandas as pd
import numpy as np

import utils

parser = argparse.ArgumentParser(description="Merge columns of multiple experiments by gene id.")

parser.add_argument('-f', '--input_file', required=True, help="Txt file with paths to expression files.")

args = parser.parse_args()


def merge_expressions(input_file):
    """Merge expressioms."""
    df = pd.DataFrame(np.nan, index=[], columns=[])
    headers = []

    with open(input_file, 'r') as tfile:
        for line in tfile:
            line = line.rstrip()
            header = os.path.basename(line).rstrip()
            if header.endswith('.tab.gz'):
                header = header[:-7]
            headers.append(header)

            with utils.gzopen(line) as csvfile:
                reader = pd.read_csv(csvfile, index_col='Gene', delimiter='\t', dtype=str)
                df = pd.concat([df, reader], axis=1)
    df.columns = headers
    return df


def write_output(input_file, df):
    """Compose a name and write output file."""
    with open(input_file, 'r') as tfile:
        first_line = next(tfile).strip()
        # Line is a path of the following form:
        # species/build/sample_slug/filetype/filename
        species, build, _, _, _ = os.path.normpath(first_line).split(os.path.sep)
        # the form of input file name is: species_build_expression-type_expressions_temp.txt
        # get Expression type from the name of input file
        new_end_name = '_'.join([input_file.split('_')[-3], 'all_expressions.tab'])
        name = os.path.join(species, build, new_end_name)

    df.to_csv(name, sep='\t', na_rep='NA', index_label='Gene')


def main():
    """Invoke when run directly as a program."""
    df = merge_expressions(args.input_file)
    write_output(args.input_file, df)


if __name__ == "__main__":
    main()

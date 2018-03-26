#!/usr/bin/env python3
"""Merge columns of multiple experiments by gene id with sample name for header."""

import argparse
import csv
import os

import pandas as pd
import numpy as np

from collections import Counter

import utils

parser = argparse.ArgumentParser(description="Merge columns of multiple experiments by gene id.")

parser.add_argument('-f', '--file_path', required=True, nargs='+', help="List with paths to expression files.")
parser.add_argument('-names', '--sample_names', required=True, nargs='+', help="List of sample names.")
parser.add_argument('-b', '--build', required=True, nargs='+', help="List of sample builds.")
parser.add_argument('-s', '--species', required=True, nargs='+', help="List of species of samples.")
parser.add_argument('-e', '--exp_type', required=True, nargs='+', help="List of expression types for samples.")


if __name__ == "__main__":

    args = parser.parse_args()

    data = {}

    for (exp_file, build, species, exp_type, sample_name) in zip(
            args.file_path, args.build, args.species, args.exp_type, args.sample_names):
        data.setdefault((build, species, exp_type), []).append([exp_file, sample_name])

    for (build, species, exp_type), data_values in data.items():
        df = pd.DataFrame(np.nan, index=[], columns=[])
        header = []
        for (exp_file, sample_name) in data_values:
            header.append(sample_name)
            with utils.gzopen(exp_file) as csvfile:
                reader = pd.read_csv(csvfile, index_col='Gene', delimiter='\t', dtype=str)
                df = pd.concat([df, reader], axis=1)

        # Add numbers to duplicated sample names.
        counts = Counter(header)
        for sample_name, num in counts.items():
            if num > 1:
                for suffix in range(1, num + 1):
                    header[header.index(sample_name)] = '{}_{}'.format(sample_name, suffix)
        df.columns = header
        name = '_'.join([species, build, exp_type, 'all_expressions.txt'])
        df.to_csv(name, sep='\t', na_rep='NA', index_label='Gene')

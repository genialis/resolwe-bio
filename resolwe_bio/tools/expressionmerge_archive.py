#!/usr/bin/env python3
"""Merge columns of multiple experiments by gene id with sample name for header."""

import argparse
from collections import Counter
import csv
import os

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Merge columns of multiple experiments by gene id.")
parser.add_argument('-f', '--file-paths', required=True, nargs='+', help="List of paths to expression files.")
parser.add_argument('-n', '--sample-names', required=True, nargs='+', help="List of sample names.")
parser.add_argument('-b', '--builds', required=True, nargs='+', help="List of sample builds.")
parser.add_argument('-s', '--species', required=True, nargs='+', help="List of species of samples.")
parser.add_argument('-e', '--exp-types', required=False, nargs='+', help="List of expression types of samples.")

na_rep = 'NA'

if __name__ == "__main__":
    args = parser.parse_args()

    # Use file in exp output field of expression data objects.
    # This is a two-column file with feature ids and expression values.
    if args.exp_types:
        data = {}

        for (exp_file, build, species, exp_type, sample_name) in zip(
                args.file_paths, args.builds, args.species, args.exp_types, args.sample_names):
            data.setdefault((build, species, exp_type), []).append([exp_file, sample_name])

        for (build, species, exp_type), data_values in data.items():
            df = pd.DataFrame(np.nan, index=[], columns=[])
            header = []
            for (exp_file, sample_name) in data_values:
                header.append(sample_name)
                reader = pd.read_csv(exp_file, index_col='Gene', delimiter='\t', dtype=str)
                df = pd.concat([df, reader], axis=1)

            # Add numbers to duplicated sample names.
            counts = Counter(header)
            for sample_name, num in counts.items():
                if num > 1:
                    for suffix in range(1, num + 1):
                        header[header.index(sample_name)] = '{}_{}'.format(sample_name, suffix)
            df.columns = header
            name = '_'.join([species, build, exp_type, 'all_expressions.txt'])
            df.to_csv(name, sep='\t', na_rep=na_rep)

    # Use file in exp_set output field of expression data objects.
    # This is a multi-column file with feature ids, gene symbols and one or more expression values.
    else:
        multi_index = ['FEATURE_ID', 'GENE_SYMBOL']
        expressions = [
            pd.read_csv(file_path, index_col=multi_index, delimiter='\t', dtype=str)
            for file_path in args.file_paths
        ]
        column_labels = {}
        for (expression, build, species) in zip(expressions, args.builds, args.species):
            for column_label in expression.columns.values:
                if column_label in multi_index:
                    continue
                column_labels.setdefault((build, species), set()).add(column_label)
        for (build, species), columns in column_labels.items():
            for column_label in columns:
                df = pd.DataFrame()
                header = []
                for i in range(len(args.file_paths)):
                    if not (args.builds[i] == build and args.species[i] == species):
                        continue
                    header.append(args.sample_names[i])
                    if column_label in expressions[i]:
                        df = pd.concat([df, expressions[i][column_label]], axis=1).drop_duplicates()
                    else:
                        empty_df = pd.DataFrame({column_label: na_rep}, index=df.index)
                        df = pd.concat([df, empty_df], axis=1).drop_duplicates()
                counts = Counter(header)
                for sample_name, num in counts.items():
                    if num > 1:
                        for suffix in range(1, num + 1):
                            header[header.index(sample_name)] = '{}_{}'.format(sample_name, suffix)
                df.columns = header
                name = '_'.join([species, build, column_label, 'all_expressions.txt'])
                df.to_csv(name, sep='\t', na_rep=na_rep)

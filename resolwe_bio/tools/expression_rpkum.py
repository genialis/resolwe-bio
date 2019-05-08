#!/usr/bin/env python3
"""Normalize expressions to RPKUM."""
import argparse
import gzip
import math
import sys

import pandas as pd

from os.path import basename
from resolwe_runtime_utils import error


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Normalize expressions using RPKUM.")
    parser.add_argument('-c', '--counts', required=True, help="Raw counts file.")
    parser.add_argument('-m', '--mappability', required=True, help="Mappability file.")
    parser.add_argument('-o', '--output', required=True, help="Output file name.")
    return parser.parse_args()


def parse_expression_file(exp_file):
    """Parse expression file to a Pandas dataframe."""
    try:
        df = pd.read_csv(
            exp_file,
            sep='\t',
            compression='gzip',
            usecols=['Gene', 'Expression'],
            index_col='Gene',
            dtype={
                'Gene': str,
                'Expression': float,
            },
        )
        return df.dropna()
    except (ValueError, OSError) as parse_error:
        print(error("Failed to read input file {}. {}".format(basename(exp_file), parse_error)))
        sys.exit(1)


def parse_mapability_file(mapability_file):
    """Parse mapability file to a Pandas dataframe."""
    try:
        mappability = pd.read_csv(
            mapability_file,
            sep='\t',
            usecols=['gene_id', 'coverage'],
            index_col='gene_id',
            dtype={
                'gene_id': str,
                'coverage': float,
            },
        )
        return mappability.dropna()
    except (ValueError, OSError) as parse_error:
        print(error("Failed to read mappability file {}. {}".format(basename(exp_file), parse_error)))
        sys.exit(1)


def get_coverage(gene_id, mappability):
    """Return gene coverage from the mapability file."""
    try:
        return mappability.loc[gene_id]['coverage']
    except KeyError:
        print(error("Feature ID {} is not present in the mappability file. "
                    "Make sure that the expressions and mappability file are "
                    "derived from the same annotations (GTF/GFF) file.".format(gene_id)))
        sys.exit(1)


def rpkum(rc, cov, lib_size):
    """Compute RPKUM."""
    if cov == 0:
        return 0
    else:
        return (math.pow(10, 9) * rc) / (lib_size * cov)


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    mappability = parse_mapability_file(args.mappability)
    expression = parse_expression_file(args.counts)
    lib_size = expression['Expression'].sum()

    # Compute RPKUM expression
    expression['Expression'] = expression.apply(
        lambda row: rpkum(row.Expression, get_coverage(row.name, mappability), lib_size),
        axis=1
    )
    expression.to_csv(args.output, columns=['Expression'], header=True, index=True, sep='\t', compression='gzip')


if __name__ == "__main__":
    main()

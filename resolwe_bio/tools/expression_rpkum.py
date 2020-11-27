#!/usr/bin/env python3
"""Normalize expressions to RPKUM."""
import argparse
import sys
from os.path import basename

import pandas as pd
from resolwe_runtime_utils import error, send_message


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Normalize expressions using RPKUM.")
    parser.add_argument("-c", "--counts", required=True, help="Raw counts file.")
    parser.add_argument("-m", "--mappability", required=True, help="Mappability file.")
    parser.add_argument("-o", "--output", required=True, help="Output file name.")
    return parser.parse_args()


def parse_expression_file(exp_file):
    """Parse expression file to a Pandas Series."""
    try:
        expression = pd.read_csv(
            exp_file,
            sep="\t",
            compression="gzip",
            usecols=["Gene", "Expression"],
            index_col="Gene",
            dtype={
                "Gene": str,
                "Expression": float,
            },
            squeeze=True,
        )
        return expression.dropna()
    except (ValueError, OSError) as parse_error:
        send_message(
            error(
                "Failed to read input file {}. {}".format(
                    basename(exp_file), parse_error
                )
            )
        )
        sys.exit(1)


def parse_mapability_file(mapability_file):
    """Parse mapability file to a Pandas Series."""
    try:
        mappability = pd.read_csv(
            mapability_file,
            sep="\t",
            usecols=["gene_id", "coverage"],
            index_col="gene_id",
            dtype={
                "gene_id": str,
                "coverage": float,
            },
            squeeze=True,
        )
        return mappability.dropna()
    except (ValueError, OSError) as parse_error:
        send_message(
            error(
                "Failed to read mappability file {}. {}".format(
                    basename(mapability_file), parse_error
                )
            )
        )
        sys.exit(1)


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    mappability = parse_mapability_file(args.mappability)
    expression = parse_expression_file(args.counts)

    missing_genes = expression.index.difference(mappability.index)
    if len(missing_genes) > 0:
        send_message(
            error(
                "Feature ID {} is not present in the mappability file. "
                "Make sure that the expressions and mappability file are "
                "derived from the same annotations (GTF/GFF) file.".format(
                    missing_genes[0]
                )
            )
        )
        sys.exit(1)

    lib_size = expression.sum()

    result = 10 ** 9 * expression / lib_size / mappability
    result[mappability == 0] = 0.0

    result.loc[expression.index].to_csv(
        args.output,
        index_label="Gene",
        header=["Expression"],
        sep="\t",
        compression="gzip",
    )


if __name__ == "__main__":
    main()

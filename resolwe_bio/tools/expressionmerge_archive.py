#!/usr/bin/env python3
"""Merge columns of multiple experiments by gene id with sample name for header."""

import argparse
import multiprocessing

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
    description="Merge columns of multiple experiments by gene id."
)
parser.add_argument(
    "-f",
    "--file-paths",
    required=True,
    nargs="+",
    help="List of paths to expression files.",
)
parser.add_argument(
    "-n", "--sample-names", required=True, nargs="+", help="List of sample names."
)
parser.add_argument(
    "-b", "--builds", required=True, nargs="+", help="List of sample builds."
)
parser.add_argument(
    "-s", "--species", required=True, nargs="+", help="List of species of samples."
)
parser.add_argument(
    "-e",
    "--exp-types",
    required=False,
    nargs="+",
    help="List of expression types of samples.",
)

na_rep = "NA"

if __name__ == "__main__":
    args = parser.parse_args()

    # Use file in exp output field of expression data objects.
    # This is a two-column file with feature ids and expression values.
    if args.exp_types:
        data = {}

        for (exp_file, build, species, exp_type, sample_name) in zip(
            args.file_paths,
            args.builds,
            args.species,
            args.exp_types,
            args.sample_names,
        ):
            data.setdefault((build, species, exp_type), []).append(
                [exp_file, sample_name]
            )

        for (build, species, exp_type), data_values in data.items():
            df = pd.DataFrame(np.nan, index=[], columns=[])
            header = []
            for (exp_file, sample_name) in data_values:
                header.append(sample_name)
                reader = pd.read_csv(
                    exp_file, index_col="Gene", delimiter="\t", dtype=str
                )
                df = pd.concat([df, reader], axis=1)

            df.columns = header
            name = "_".join([species, build, exp_type, "all_expressions.txt"])
            df.to_csv(name, sep="\t", na_rep=na_rep)

    # Use file in exp_set output field of expression data objects.
    # This is a multi-column file with feature ids, gene symbols and one or more expression values.
    else:
        multi_index = ["FEATURE_ID", "GENE_SYMBOL"]

        def read_csv(file_path):
            """Return DataFrame representation of CSV file."""
            return pd.read_csv(
                file_path, index_col=multi_index, delimiter="\t", dtype=str
            )

        expressions = multiprocessing.Pool().map(read_csv, args.file_paths)

        items = []
        for (species, build, expression) in zip(args.species, args.builds, expressions):
            for column_label in expression.columns:
                if column_label in multi_index:
                    continue
                item = {
                    "species": species,
                    "build": build,
                    "column_label": column_label,
                }
                if item not in items:
                    items.append(item)

        def merge_expressions(item):
            """Merge expressions of multiple samples and save them into CSV file."""
            dfs = [
                expressions[i][item["column_label"]]
                for i in range(len(args.file_paths))
                if (
                    args.species[i] == item["species"]
                    and args.builds[i] == item["build"]
                    and item["column_label"] in expressions[i].columns
                )
            ]
            df = pd.concat(dfs, axis=1)

            header = [
                args.sample_names[i]
                for i in range(len(args.file_paths))
                if (
                    args.species[i] == item["species"]
                    and args.builds[i] == item["build"]
                    and item["column_label"] in expressions[i].columns
                )
            ]
            df.columns = header

            name = "_".join(
                [
                    item["species"],
                    item["build"],
                    item["column_label"],
                    "all_expressions.txt",
                ]
            )
            df.to_csv(name, sep="\t", na_rep=na_rep)

        multiprocessing.Pool().map(merge_expressions, items)

#!/usr/bin/env python3
"""Parse Diff Exp output files."""

import argparse
import json

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from resolwe_runtime_utils import error


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse Diff Exp output files")
    parser.add_argument("raw_file", help="DE analysis output file (.tab).")
    parser.add_argument("output_json", help="Output JSON")
    parser.add_argument("output_file", help="Output file")
    parser.add_argument("--gene_id", help="Gene_IDs column name", type=str)
    parser.add_argument("--fdr", help="FDR column name", type=str)
    parser.add_argument("--pvalue", help="Pvalue column name", type=str)
    parser.add_argument("--fwer", help="FWER column name", type=str)
    parser.add_argument("--logodds", help="Log Odds column name", type=str)
    parser.add_argument("--logfc", help="logfc column name", type=str)
    parser.add_argument("--stat", help="Statistics column name", type=str)
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    de_data = pd.read_csv(args.raw_file, sep="\t")
    de_data.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
    de_data.fillna(value=1, inplace=True)
    columns = {}
    col_order = []

    # Make sure all listed numeric columns are valid numeric variables based
    # on a union of numeric column names from cuffdiff, edgeR, deseq2 and test
    # files.
    numeric_columns = [
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "value_1",
        "value_2",
        "log2(fold_change)",
        "test_stat",
        "p_value",
        "q_value",
        "logfc",
        "fdr",
        "stat",
        "logFC",
        "logCPM",
        "LR",
        "Pvalue",
        "FDR",
    ]
    de_columns = de_data.columns

    for column in numeric_columns:
        if column not in de_columns:
            continue

        if not is_numeric_dtype(de_data[column]):
            msg = (
                f"Column {column} is not numeric. Please make sure "
                f"that the input file has valid numeric values (i.e. "
                f"periods for decimal places)."
            )
            print(error(msg))
            raise ValueError(msg)

    if args.gene_id:
        if args.gene_id == "index":
            columns["gene_id"] = list(de_data.index.astype(str))
            col_order.append("gene_id")
        else:
            columns["gene_id"] = list(de_data[args.gene_id].astype(str))
            col_order.append("gene_id")

    if args.logfc:
        col = np.array(de_data[args.logfc])
        col[np.isinf(col)] = 0
        columns["logfc"] = list(col)
        col_order.append("logfc")

    if args.fdr:
        columns["fdr"] = list(de_data[args.fdr])
        col_order.append("fdr")

    if args.pvalue:
        columns["pvalue"] = list(de_data[args.pvalue])
        col_order.append("pvalue")

    if args.fwer:
        columns["fwer"] = list(de_data[args.fwer])
        col_order.append("fwer")

    if args.logodds:
        columns["logodds"] = list(de_data[args.logodds])
        col_order.append("logodds")

    if args.stat:
        columns["stat"] = list(de_data[args.stat])
        col_order.append("stat")

    with open(args.output_json, "w") as f:
        json.dump(columns, f, separators=(",", ":"), allow_nan=False)

    outdf = pd.DataFrame(columns)
    outdf = outdf[col_order]
    outdf.to_csv(args.output_file, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()

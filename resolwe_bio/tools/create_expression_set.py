#!/usr/bin/env python3
"""Map Feature IDs to Feature Symbols and write to the expression set file."""

import argparse
import gzip
import json
import sys

import pandas as pd
import resdk
from resolwe_runtime_utils import error, send_message, warning

CHUNK_SIZE = 1000


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Append column with feature symbols to expressions file."
    )
    parser.add_argument(
        "--expressions", required=True, help="Expression file with Feature IDs."
    )
    parser.add_argument(
        "--source_db", type=str, required=True, help="Feature IDs source database."
    )
    parser.add_argument(
        "--species", type=str, required=True, help="Feature IDs species."
    )
    parser.add_argument(
        "--output_name", type=str, required=True, help="Output file name prefix."
    )
    parser.add_argument(
        "--expressions_type", type=str, default="RAW_COUNT", help="Type of expression."
    )
    parser.add_argument(
        "--norm_expressions",
        nargs="+",
        help="Additional (normalized) expression files.",
    )
    parser.add_argument(
        "--norm_expressions_type",
        nargs="+",
        type=str,
        help="Type of additional expression files.",
    )
    return parser.parse_args()


def parse_expression_file(exp_file, exp_type):
    """Parse expression file to a Pandas dataframe."""
    with gzip.open(exp_file) as exp:
        df = pd.read_csv(exp, sep="\t")

        ALLOWED_COLUMNS = ["Gene", "Transcript", "Expression"]
        if not all(
            column_label in ALLOWED_COLUMNS for column_label in df.columns.values
        ):
            send_message(
                error(
                    "Invalid column headers {} in file {}.".format(
                        df.columns.values, exp_file
                    )
                )
            )
            sys.exit(1)

        df.rename(
            index=str,
            columns={
                "Gene": "FEATURE_ID",
                "Transcript": "FEATURE_ID",
                "Expression": exp_type,
            },
            inplace=True,
        )
        # Cast FEATURE_ID column to string
        df["FEATURE_ID"] = df["FEATURE_ID"].astype("str")
        # Remove any possible empty rows from the input file
        df.dropna(inplace=True)

    return df


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    if args.norm_expressions and args.norm_expressions_type:
        if len(args.norm_expressions) != len(args.norm_expressions_type):
            send_message(
                error(
                    "The number of additional expression files must match the number of specified "
                    "expressions types."
                )
            )
            sys.exit(1)

    if args.norm_expressions_type:
        exp_types = [args.expressions_type] + args.norm_expressions_type
        if len(exp_types) != len(set(exp_types)):
            send_message(
                error(
                    "The union of the main expression type ({}) and additional normalized expression types {} "
                    "does not contain unique items.".format(
                        args.expressions_type, args.norm_expressions_type
                    )
                )
            )
            sys.exit(1)

    res = resdk.Resolwe()

    feature_dict = {}
    df = parse_expression_file(args.expressions, args.expressions_type)

    # Get a list of feature IDs
    input_features = df["FEATURE_ID"].tolist()

    # Split feature IDs into chunks with max size of 10000 elements
    features_sublists = [
        input_features[i : i + CHUNK_SIZE]
        for i in range(0, len(input_features), CHUNK_SIZE)
    ]

    # Fetch features from KB and add them to {feature_id: feature_name} mapping dict
    for fsublist in features_sublists:
        features = res.feature.filter(
            source=args.source_db, species=args.species, feature_id__in=fsublist
        )
        feature_dict.update({f.feature_id: f.name for f in features})

    # Map gene symbols to feature IDs
    df["GENE_SYMBOL"] = df["FEATURE_ID"].map(feature_dict)

    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        send_message(
            warning(
                "{} feature(s) could not be mapped to the associated feature symbols.".format(
                    sum(df.isnull().values.ravel())
                )
            )
        )

    # Merge additional expression files with the original data frame
    if args.norm_expressions and args.norm_expressions_type:
        for exp_file, exp_type in zip(
            args.norm_expressions, args.norm_expressions_type
        ):
            exp_df = parse_expression_file(exp_file, exp_type)
            df = df.merge(exp_df, on="FEATURE_ID")

    # Reorder the columns in dataframe
    columns = ["FEATURE_ID", "GENE_SYMBOL", args.expressions_type]
    if args.norm_expressions_type:
        columns = columns + args.norm_expressions_type
    df = df[columns]

    # Replace NaN values with empty string
    df.fillna("", inplace=True)

    # Write to file
    df.to_csv(
        args.output_name + ".txt.gz",
        header=True,
        index=False,
        sep="\t",
        compression="gzip",
    )

    # Write to JSON
    df_dict = df.set_index("FEATURE_ID").to_dict(orient="index")
    with open(args.output_name + ".json", "w") as f:
        json.dump({"genes": df_dict}, f, allow_nan=False)


if __name__ == "__main__":
    main()

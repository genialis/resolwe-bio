#!/usr/bin/env python3
"""Validate the uploaded master file."""

import argparse
import csv

from resolwe_runtime_utils import error, send_message


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run gene ontology enrichment analysis."
    )
    parser.add_argument("master_file", type=str, help="Master file.")
    return parser.parse_args()


def check_dna_sequence(sequence):
    """Check if a given sequence contains only the allowed letters A, C, T, G."""
    return len(sequence) != 0 and all(
        base.upper() in ["A", "C", "T", "G"] for base in sequence
    )


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    amplicon_names = set()

    with open(args.master_file, newline="") as masterfile:
        reader = csv.reader(masterfile, delimiter="\t")
        for row in reader:
            if len(row) != 12:
                send_message(
                    error("Uploaded master file must contain exactly 12 columns.")
                )
            if not check_dna_sequence(row[10]):
                send_message(error("11th column must contain a DNA sequence."))
            if not check_dna_sequence(row[11]):
                send_message(error("12th column must contain a DNA sequence."))

            amp_name = row[3]
            if amp_name not in amplicon_names:
                amplicon_names.add(amp_name)
            else:
                send_message(
                    error(
                        "Amplicon names must be unique. Amplicon {} is seen multiple times.".format(
                            amp_name
                        )
                    )
                )


if __name__ == "__main__":
    main()

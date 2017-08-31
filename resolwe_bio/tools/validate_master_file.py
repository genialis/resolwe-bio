#!/usr/bin/env python3
"""Validate the uploaded master file."""

import argparse
import csv


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run gene ontology enrichment analysis.")
    parser.add_argument('master_file', type=str, help="Master file.")
    return parser.parse_args()


def check_dna_sequence(sequence):
    """Check if a given sequence contains only the allowed letters A, C, T, G."""
    return len(sequence) != 0 and all(base in ['A', 'C', 'T', 'G'] for base in sequence)


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    with open(args.master_file, newline='') as masterfile:
        reader = csv.reader(masterfile, delimiter='\t')
        for row in reader:
            if len(row) != 12:
                print('{"proc.error":"Uploaded master file must contain exactly 12 columns."}')
                exit(1)
            if not check_dna_sequence(row[10]):
                print('{"proc.error":"11th column must contain a DNA sequence."}')
                exit(1)
            if not check_dna_sequence(row[11]):
                print('{"proc.error":"12th column must contain a DNA sequence."}')
                exit(1)


if __name__ == "__main__":
    main()

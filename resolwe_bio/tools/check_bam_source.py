#!/usr/bin/python3
"""Check if the source of bam file were rRNA or globin RNA."""

import argparse
import os
import sys

GLOBIN_RRNA_DIR = "/opt/chrom_mappings/globin_rrna/"
# These are files that contain names of globin RNA and rRNA sequences."
GLOBIN_RRNA_FILES = {
    "Homo sapiens": [
        "Homo_sapiens_rRNA.txt",
        "Homo_sapiens_globin_RNA.txt",
    ],
    "Mus musculus": [
        "Mus_musculus_rRNA.txt",
        "Mus_musculus_globin_RNA.txt",
    ],
    "Rattus norvegicus": [
        "Rattus_norvegicus_rRNA.txt",
        "Rattus_norvegicus_globin_RNA.txt",
    ],
}


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", help="File with chromosome names from bam file.")
    parser.add_argument("--species", help="Species")
    return parser.parse_args()


def parse_chrom_names(filename):
    """Parse file with chromosome names."""
    chrom_names = set()
    with open(filename) as handle:
        for line in handle:
            line = line.strip().split("\t")
            chrom_names.add(line[0])
    return chrom_names


def get_globin_rrna_chrom_names(species):
    """Get chromosome names of globin RNA and rRNA."""
    if species not in GLOBIN_RRNA_FILES:
        sys.exit(0)

    chrom_names = set()
    for basename in GLOBIN_RRNA_FILES[species]:
        filename = os.path.join(GLOBIN_RRNA_DIR, basename)
        chrom_names.update(parse_chrom_names(filename))
    return chrom_names


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    globin_rrna_chrom_names = get_globin_rrna_chrom_names(
        args.species.strip("'").capitalize()
    )
    bam_chrom_names = parse_chrom_names(args.infile)

    if bam_chrom_names.issubset(globin_rrna_chrom_names):
        print("globin_or_rrna")


if __name__ == "__main__":
    main()

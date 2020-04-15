#!/usr/bin/env python3
"""Convert ggf3 to gtf file."""
from __future__ import absolute_import, division, print_function

import argparse


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Gff3 to gtf convertion")
    parser.add_argument("gff3_file", help="gff3 file")
    parser.add_argument("gtf_file", help="output file (gtf)")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    inFile = open(args.gff3_file, "r")

    with open(args.gtf_file, "w") as f:
        for line in inFile:
            line = line.strip()
            if line[0] != "#":
                # split line into columns by tab
                data = line.split("\t")

                ID = ""
                if data[2] == "gene":  # if the feature is a gene
                    # get the id
                    ID = data[-1].split("ID=")[-1].split(";")[0]

                else:  # if the feature is anything else
                    # get the parent as the ID
                    ID = data[-1].split("Parent=")[-1].split(";")[0]

                # modify the last column
                data[-1] = 'gene_id "{}"; transcript_id "{}";'.format(ID, ID)

                # print out this new GTF line
                gtf = "\t".join(data) + "\n"

                f.write(gtf)


if __name__ == "__main__":
    main()

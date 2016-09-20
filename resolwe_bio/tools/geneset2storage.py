#!/usr/bin/env python2
"""Parse Gene set file."""
from __future__ import absolute_import, division, print_function

import argparse
import json


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse Gene set file")
    parser.add_argument('geneset_file', help="Gene set file (.tab).")
    parser.add_argument('output_file', help='Output JSON file')
    return parser.parse_args()


def main():
    """Invoked when run directly as a program."""
    args = parse_arguments()

    with open(args.geneset_file) as infile:
        geneset = [gene.strip() for gene in infile]

        with open(args.output_file, 'w') as f:
            json.dump({'genes': geneset}, f, separators=(',', ':'), allow_nan=False)


if __name__ == "__main__":
    main()

#!/usr/bin/env python2
"""Parse Gene set file."""

import argparse
import json


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse Gene set file")
    parser.add_argument('geneset_file', help="Gene set file (.tab).")
    return parser.parse_args()


def main():
    """Invoked when run directly as a program."""
    args = parse_arguments()

    with open(args.geneset_file) as infile:
        geneset = [gene.strip() for gene in infile]
        data = {'geneset_json': {'geneset_json': geneset}}
        print(json.dumps(data, separators=(',', ':'), allow_nan=False))


if __name__ == "__main__":
    main()

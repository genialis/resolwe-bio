#!/usr/bin/env python3
"""Parse Gene set file."""

import argparse
import gzip
import json


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse Gene set file")
    parser.add_argument("geneset_file", help="Gene set file (.tab).")
    parser.add_argument("output_json", help="Output JSON file")
    parser.add_argument("output_file", help="Output file")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    with open(args.geneset_file, "rU") as infile:
        # skip empty lines in input gene set file
        genes = [str(line.strip()) for line in infile if line.strip()]
        geneset = sorted(set(genes))

        if len(genes) != len(geneset):
            print('{"proc.warning":"Removed duplicated genes."}')

        with open(args.output_json, "w") as json_out:
            json.dump(
                {"genes": geneset}, json_out, separators=(",", ":"), allow_nan=False
            )

        with gzip.open(args.output_file, "w") as file_out:
            file_out.write("\n".join(geneset).encode("utf-8"))


if __name__ == "__main__":
    main()

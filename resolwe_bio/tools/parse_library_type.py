#!/usr/bin/env python3
"""Return library type information from the Salmon output file."""
import argparse
import json

from resolwe_runtime_utils import error, save, send_message


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse library type information.")
    parser.add_argument("input_file", help="Salmon library type information file.")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    with open(args.input_file) as infile:
        data = json.load(infile)
        if "expected_format" in data and "compatible_fragment_ratio" in data:
            send_message(save("strandedness", data["expected_format"]))
            send_message(
                save("fragment_ratio", str(round(data["compatible_fragment_ratio"], 2)))
            )
        else:
            send_message(error("Cannot parse library type information file."))


if __name__ == "__main__":
    main()

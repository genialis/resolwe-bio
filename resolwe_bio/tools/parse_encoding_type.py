#!/usr/bin/env python3
"""Parse encoding type."""
import argparse

parser = argparse.ArgumentParser(description="Parse encoding type.")

parser.add_argument("report_file", help="FastQC report file")
args = parser.parse_args()

encoding = ""

with open(args.report_file) as report:
    for line in report:
        if line.startswith("Encoding"):
            encoding = line.strip().split("\t")[1]
            break
    if encoding != "":
        print(encoding)
    else:
        print("Unknown")

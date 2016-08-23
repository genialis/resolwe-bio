#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Add header row."""
import csv
import argparse

parser = argparse.ArgumentParser(description="Add header row")
parser.add_argument('in_file', help="In-file")
parser.add_argument('out_file', help="Out-file")
args = parser.parse_args()

with open(args.out_file, 'wb') as outcsv:
    writer = csv.writer(outcsv, delimiter='\t')
    writer.writerow(["Gene", "Expression"])

    with open(args.in_file, "r+") as incsv:
        reader = csv.reader(incsv, delimiter="\t")
        writer.writerows(row for row in reader)

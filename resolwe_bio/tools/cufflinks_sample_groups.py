#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Prepare groups and replicates info."""
from __future__ import absolute_import, division, print_function

import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="Prepare groups and replicates info")
parser.add_argument('samples', help="All samples")
parser.add_argument('replicates', help="Define replicates and groups")
args = parser.parse_args()

samples = args.samples.split(",")
replicates = args.replicates.split(",")

d = defaultdict(list)

# This sort must produce the same order as cuffnorm_output_table_headers.py
for s, r in zip(samples, replicates):
    d[r].append(s)

return_string = ""

for x in sorted(d):
    return_string = return_string + ",".join(d[x]) + " "

print(return_string)

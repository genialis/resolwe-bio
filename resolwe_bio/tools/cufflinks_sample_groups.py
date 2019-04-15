#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Prepare groups and replicates info."""

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Prepare groups and replicates info")
parser.add_argument('samples', help="All samples")
parser.add_argument('replicates', help="Define replicates and groups")
args = parser.parse_args()

samples = args.samples.split(",")
replicates = args.replicates.split(",")

d = defaultdict(list)
for s, r in zip(samples, replicates):
    d[r].append(s)

return_string = ""
for x in sorted(d):
    sample_groups = ','.join(d[x])
    return_string = return_string + "'{}' ".format(sample_groups)

print(return_string)

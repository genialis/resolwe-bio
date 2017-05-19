#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Create group labels from replicates."""
from __future__ import absolute_import, division, print_function

import argparse

parser = argparse.ArgumentParser(description="Create group labels")
parser.add_argument('replicates', help="Define replicates and groups")
args = parser.parse_args()

replicates = set(args.replicates.split(","))

labels = ""
for rep in replicates:
    labels = labels + "Group" + rep + ","

print(labels)

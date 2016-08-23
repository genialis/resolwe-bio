#!/usr/bin/env python
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Orthologues file upload verification."""
from __future__ import absolute_import, division, print_function

import argparse
import os


parser = argparse.ArgumentParser(description="Orthologues file upload verification")
parser.add_argument('orthologues', help="Orthologues file")
args = parser.parse_args()

if not os.path.isfile(args.orthologues):
    raise ValueError("File {} does not exist".format(args.orthologues))

orthologues = open(args.orthologues)

for line in orthologues:
    if not len(line.strip().split("\t")) == 2 or line.startswith("#"):
        print("{'proc.rc':1}")
        exit(1)

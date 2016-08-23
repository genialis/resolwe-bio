#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Check if sample names are unique."""
from __future__ import absolute_import, division, print_function

import argparse
import json


parser = argparse.ArgumentParser(description="Check if sample names are unique")
parser.add_argument('samples', help="All samples")
args = parser.parse_args()

samples = args.samples.split(",")

if len(samples) > len(set(samples)):
    print(json.dumps({"proc.error": "Sample names must be unique."}))

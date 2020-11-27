#!/usr/bin/env python3
"""Check if sample names are unique."""
import argparse

from resolwe_runtime_utils import error, send_message

parser = argparse.ArgumentParser(description="Check if sample names are unique")
parser.add_argument("samples", help="All samples")
args = parser.parse_args()

samples = args.samples.split(",")

if len(samples) > len(set(samples)):
    send_message((error("Sample names must be unique.")))

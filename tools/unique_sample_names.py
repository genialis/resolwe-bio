#!/usr/bin/env python2
import argparse
import json

parser = argparse.ArgumentParser(description="Check if sample names are unique")
parser.add_argument('samples', help="All samples")
args = parser.parse_args()

samples = args.samples.split(",")

if len(samples) > len(set(samples)):
    print json.dumps({"proc.error": "Sample names must be unique."})

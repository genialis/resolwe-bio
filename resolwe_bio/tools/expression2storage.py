#!/usr/bin/env python3
# XXX: Refactor to a comand line tool and remove pylint disable
"""Save expressions to storage."""

import argparse
import json
import os

import utils

parser = argparse.ArgumentParser(description="Parses expressions for storage.")
parser.add_argument("input", help="Input expression file")
parser.add_argument("--output", help="Output JSON file")
args = parser.parse_args()

if not args.input:
    print('{"rc":"1"}')
    exit(1)

if not (args.input and os.path.isfile(args.input)):
    print('{"rc":"1"}')
    exit(1)


def isfloat(value):
    """Check if value is float."""
    try:
        float(value)
        return True
    except ValueError:
        return False


with utils.gzopen(args.input) as f:
    # Split lines by tabs
    # Ignore lines without a number in second column
    # Build a dictionary of gene-expression pairs
    exp = {
        "genes": {
            gene_exp[0]: float(gene_exp[1])
            for gene_exp in (l.split("\t") for l in f)
            if len(gene_exp) == 2 and isfloat(gene_exp[1])
        }
    }

if args.output:
    with open(args.output, "w") as f:
        json.dump(exp, f)
else:
    print('{"exp_json":%s}' % json.dumps(exp, separators=(",", ":")))

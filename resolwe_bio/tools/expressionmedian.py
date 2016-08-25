#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Median gene expressions of multiple experiments."""
from __future__ import absolute_import
import argparse
import collections
import csv
import os
import sys

import numpy as np  # pylint: disable=import-error

import utils


parser = argparse.ArgumentParser(description='Median gene expressions of multiple experiments.')
parser.add_argument('files', nargs='*', help='expression files')
parser.add_argument('--name', help='expression column name')
parser.add_argument('--out', help='output file')

args = parser.parse_args()

expressions = collections.defaultdict(list)
for f in args.files:
    if not os.path.isfile(f):
        exit(1)

    base, ext = os.path.splitext(f)
    delimiter = ';' if ext == '.csv' else '\t'

    with utils.gzopen(f) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        header = reader.next()
        for gene, exp in reader:
            expressions[gene].append(float(exp))

genes = sorted(expressions.keys())
medians = [np.median(expressions[g]) for g in genes]

fhandler = open(args.out, 'w') if args.out else sys.stdout

fhandler.write('Gene\t{}Median\n'.format(args.name if args.name else ''))
for gene, med in zip(genes, medians):
    fhandler.write('{}\t{:.6f}\n'.format(gene, med))

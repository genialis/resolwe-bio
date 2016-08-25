#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Merge columns of multiple experiments by gene id."""
from __future__ import absolute_import
import argparse
import csv
import os
import sys

from itertools import chain

import utils


parser = argparse.ArgumentParser(description='Merge columns of multiple experiments by gene id.')
parser.add_argument('files', nargs='*', help='expression files')
parser.add_argument('--experiments', nargs='+', help='experiment ids')
parser.add_argument('--genes', nargs='+', help='filter genes')
parser.add_argument('--intersection', action='store_true', help='merge by intersection of gene ids')
parser.add_argument('--out', help='output file')


args = parser.parse_args()

# if args.experiments and len(args.experiments) != len(args.files):
#     raise ValueError("Number of experiments must match the number of files")


genes = set()
expressions = []
headers = []
op = set.intersection if args.intersection else set.union
offset = 0

for f in args.files:
    if not os.path.isfile(f):
        exit(1)

    base, ext = os.path.splitext(f)
    delimiter = ';' if ext == '.csv' else '\t'

    with utils.gzopen(f) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        header = reader.next()[1:]
        headers.append(args.experiments[offset:offset + len(header)] if args.experiments else header)
        offset += len(headers[-1])
        expressions.append(dict((r[0], r[1:]) for r in reader))
        genes = set(expressions[-1].keys()) if args.intersection and not genes else op(genes, expressions[-1].keys())

if args.genes:
    genes = genes.intersection(args.genes)

genes = sorted(genes)
he = zip(headers, expressions)
rows = [dict(chain.from_iterable([zip(h, e[g]) for h, e in he if g in e]), **{'Gene': g}) for g in genes]
fhandler = open(args.out, 'wb') if args.out else sys.stdout

writer = csv.DictWriter(fhandler, ['Gene'] + [h for subheader in headers for h in subheader], delimiter='\t')
writer.writeheader()
writer.writerows(rows)

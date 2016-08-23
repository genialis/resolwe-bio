#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Change genes names to orthologues ones."""
from __future__ import absolute_import, division, print_function

import argparse
import csv

import utils


parser = argparse.ArgumentParser(description='Change genes names to orthologues ones.')

parser.add_argument('ortholog_file', help='file with orthologues')
parser.add_argument('genes', nargs='*', help='genes names')

args = parser.parse_args()

orthologues = {}

with utils.gzopen(args.ortholog_file) as ortholog_tsv:
    for ortholog in csv.reader(ortholog_tsv, delimiter='\t'):
        orthologues[ortholog[0]] = ortholog[1]

genes = args.genes
for i, gene in enumerate(genes):
    if gene in orthologues:
        genes[i] = orthologues[gene]

print(' '.join(genes))

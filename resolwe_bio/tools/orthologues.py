#!/usr/bin/env python2
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
for i in range(len(genes)):
    if genes[i] in orthologues:
        genes[i] = orthologues[genes[i]]

print ' '.join(genes)

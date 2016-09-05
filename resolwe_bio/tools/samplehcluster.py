#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering of samples."""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import json
import csv
import gzip

import numpy as np
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage

from resolwe_runtime_utils import error, warning


parser = argparse.ArgumentParser(description='Hierarchical clustering.')

parser.add_argument('sample_files', nargs='+', help='sample files')
parser.add_argument('-s', '--sampleids', nargs='+', default=[], help='sample ids')
parser.add_argument('-g', '--genes', nargs='+', default=[], help='subset of gene ids')
parser.add_argument('-d', '--dstfunc', default='euclidean', help='distance function')
parser.add_argument('-l', '--linkage', default='average', help='clustering linkage function')
parser.add_argument('--filter', help="Filter genes with low expression", action="store_true")

args = parser.parse_args()

distance_map = {
    'spearman': lambda x, y: 1 - spearmanr(x, y).correlation,
    'pearson': lambda x, y: 1 - np.corrcoef(x, y)[0][1],
    'euclidean': 'euclidean'
}

if args.dstfunc not in distance_map:
    msg = "Invalid distance function {}".format(args.dstfunc)
    print(error(msg))
    raise ValueError(msg)

if args.linkage not in ['average', 'single', 'complete']:
    msg = "Invalid clustering linkage function {}".format(args.linkage)
    print(error(msg))
    raise ValueError(msg)

if not args.sampleids or len(args.sampleids) != len(args.sample_files):
    msg = "Number of sample ids must match the number of files"
    print(error(msg))
    raise ValueError(msg)

# read data
matrix = []
for fname in args.sample_files:
    matrix.append([])

    with gzip.open(fname) as csvfile:
        sample_reader = csv.reader(csvfile, delimiter=b'\t')
        next(sample_reader)  # Skip the header

        for row in sample_reader:
            gene_name, gene_value = row

            if not args.genes or gene_name in args.genes:
                matrix[-1].append(float(gene_value))

matrix = np.array(matrix)

matrix_sum = np.sum(matrix, axis=0)  # sum of expressions for each gene
genes_zero = np.where(matrix_sum < 0.1)[0]

if args.filter:
    matrix = np.delete(matrix, genes_zero, axis=1)

    if matrix.shape[1] == 0:
        raise ValueError("Expressions of all selected genes are 0")

distance = distance_map[args.dstfunc.lower()]
cluster = linkage(matrix, method=args.linkage.lower(), metric=distance)

distance_sum = cluster[:, 2].sum()
if distance_sum < 0.1:
    print(warning('All sample distances are 0.'))

dend = dendrogram(cluster, no_plot=True)

sample_ids = {}
for i, sample_id in enumerate(args.sampleids):
    sample_ids[i + 1] = {'id': sample_id}

output = {'cluster': {'linkage': cluster.tolist(),
                      'samples_names': sample_ids,
                      'order': dend['leaves']}}

print(json.dumps(output, separators=(',', ':')))

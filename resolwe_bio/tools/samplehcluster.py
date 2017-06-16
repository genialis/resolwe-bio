#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering of samples."""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import json
import csv
import gzip

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage

from resolwe_runtime_utils import error, warning


parser = argparse.ArgumentParser(description='Hierarchical clustering.')

parser.add_argument('sample_files', nargs='+', help='sample files')
parser.add_argument('-s', '--sampleids', nargs='+', default=[], help='sample ids')
parser.add_argument('-g', '--genes', nargs='+', default=[], help='subset of gene ids')
parser.add_argument('-d', '--dstfunc', help='distance function')
parser.add_argument('-l', '--linkage', help='clustering linkage function')
parser.add_argument('--log2', action='store_true', help='Log2')
parser.add_argument('--output', help='Output JSON file')

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
frame = pd.DataFrame()
matrix = []
for fname in args.sample_files:
    df = pd.read_csv(fname, sep='\t', header=0, index_col=0, compression='gzip')
    matrix.append(df)
frame = pd.concat(matrix, axis=1)
frame.index = frame.index.map(str)

# list of genes from files
genes_from_files = list(frame.index.values)
# list of genes defined in args.genes
gene_subset = args.genes

# insert genes that are in args.genes and are missing in all samples
for gene in gene_subset:
    if gene not in genes_from_files:
        frame.loc[gene] = [np.NaN] * len(args.sample_files)

    frame = frame.loc[gene_subset]  # select only genes from args.genes

# list of all genes that are in data frame (from files and in args.genes)
all_genes = list(frame.index.values)

matrix = frame.fillna(value=0)
matrix = np.array(matrix)
# samples are in rows
matrix = np.transpose(matrix)

# sum of  all expressions per sample to find samples without any exprerssed genes
sample_sum = np.sum(matrix, axis=1)
sample_zero = np.where(sample_sum < 0.1)[0]

# list of samples with zero expressions
list_zero_samples = [int(args.sampleids[index]) for index in sample_zero]

# sum of  all expressions per gene to find genes that are not exprerssed in any sample
genes_sum = np.sum(matrix, axis=0)
genes_zero = np.where(genes_sum < 0.1)[0]

# delete samples where none of the genes are expressed
matrix = np.delete(matrix, sample_zero, axis=0)
if matrix.shape[0] == 0:
    msg = "Expressions of selected genes are 0. Please select different samples or genes."
    print(error(msg))
    raise ValueError(msg)

if args.log2:
    matrix = np.log2(matrix + 1.0)

distance = distance_map[args.dstfunc]
cluster = linkage(matrix, method=args.linkage, metric=distance)

distance_sum = cluster[:, 2].sum()
if distance_sum < 0.1:
    msg = 'All sample distances are 0.'
    print(warning(msg))

dend = dendrogram(cluster, no_plot=True)

# create list of samples without deleted samples (with zero expression for all genes)
calculated_samples = [sample for sample in args.sampleids if sample not in list_zero_samples]

calculated_sample_ids = {}
for i, sample_id in enumerate(calculated_samples):
    calculated_sample_ids[i] = {'id': int(sample_id)}

cluster = {
    'linkage': cluster.tolist(),
    'sample_ids': calculated_sample_ids,
    'order': dend['leaves'],
    'zero_gene_symbols': [all_genes[index] for index in genes_zero],
    'zero_sample_ids': list_zero_samples
}

if args.output:
    with open(args.output, 'w') as f:
        json.dump(cluster, f)
else:
    print(json.dumps({'cluster': cluster}, separators=(',', ':')))

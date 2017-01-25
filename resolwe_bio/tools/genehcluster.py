#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering of genes."""
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
parser.add_argument('-d', '--dstfunc', default='euclidean', help='distance function')
parser.add_argument('-l', '--linkage', default='average', help='clustering linkage function')

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

# sum of expressions for each sample
sample_sum = np.sum(matrix, axis=0)
sample_zero = np.where(sample_sum < 0.1)[0]

# sum of  all expressions per gene to find genes that are not exprerssed in any sample
genes_sum = np.sum(matrix, axis=1)
genes_zero = np.where(genes_sum < 0.1)[0]
# list of genes that are not exprerssed in any sample
list_zero_genes = [all_genes[index] for index in genes_zero]

# delete genes that have zero expression in all samples
matrix = np.delete(matrix, genes_zero, axis=0)
if matrix.shape[0] == 0:
    msg = "Expressions of selected genes are 0. Please selectdifferent samples or genes."
    print(error(msg))
    raise ValueError(msg)

# calculation of hierarhical clustering
distance = distance_map[args.dstfunc]
cluster = linkage(matrix, method=args.linkage, metric=distance)

# print warning if all distances are 0
distance_sum = cluster[:, 2].sum()
if distance_sum < 0.1:
    msg = 'All genes distances are 0.'
    print(warning(msg))

dend = dendrogram(cluster, no_plot=True)

# create list of genes without deleted genes (with zero expression in all samples)
calculated_genes = [gene for gene in all_genes if gene not in list_zero_genes]

calculated_genes_dict = {}
for i, gene in enumerate(calculated_genes):
    calculated_genes_dict[i] = {'gene': (gene)}

output = {
    'cluster': {
        'linkage': cluster.tolist(),
        'gene_symbols': calculated_genes_dict,
        'order': dend['leaves'],
        'zero_gene_symbols': list_zero_genes,
        'zero_sample_ids': [int(args.sampleids[index]) for index in sample_zero]
    }
}

print(json.dumps(output, separators=(',', ':')))

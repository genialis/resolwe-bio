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
from scipy.spatial.distance import pdist

from resolwe_runtime_utils import error, warning


parser = argparse.ArgumentParser(description='Hierarchical clustering.')

parser.add_argument('sample_files', nargs='+', help='sample files')
parser.add_argument('-s', '--sampleids', nargs='+', default=[], help='sample ids')
parser.add_argument('-g', '--genes', nargs='+', default=[], help='subset of gene ids')
parser.add_argument('-d', '--dstfunc', default='euclidean', help='distance function')
parser.add_argument('-l', '--linkage', default='average', help='clustering linkage function')
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
    msg = "Expressions of selected genes are 0. Please select different samples or genes."
    print(error(msg))
    raise ValueError(msg)

if args.log2:
    matrix = np.log2(matrix + 1.0)

# calculation of hierarhical clustering
distance = distance_map[args.dstfunc]
cluster = linkage(matrix, method=args.linkage, metric=distance)

# print warning if all distances are 0
distance_sum = cluster[:, 2].sum()
if distance_sum < 0.1:
    msg = 'All genes distances are 0.'
    print(warning(msg))


def swap_children(Z, i, j):
    """Swap children."""
    if int(Z[i][j]) >= len(Z) + 1:
        row = Z[int(Z[i][j]) - len(Z) - 1]
        row[0], row[1] = row[1], row[0]
        for k in range(2):
            swap_children(Z, int(Z[i][j]) - len(Z) - 1, k)


def get_condensed_from_square(i, j, n):
    """Get condensed index from square indices."""
    assert i != j, "Error: no diagonal elements in condensed matrix."
    if i < j:
        i, j = j, i
    return int(n * j - j * (j + 1) / 2 + i - 1 - j)


def get_child(Z, i, j, k):
    """Get child."""
    if Z[i][j] < len(Z) + 1:
        return int(Z[i][j])
    return get_child(Z, int(Z[i][j]) - len(Z) - 1, k, k)


def get_children(Z, i):
    """Get children."""
    # ll lr rl rr
    return [[get_child(Z, i, j, k) for k in range(2)] for j in range(2)]


def is_ordered(Z, D, i):
    """Check if children are ordered."""
    children = get_children(Z, i)
    d = [D[get_condensed_from_square(children[0][j], children[1][k], len(Z) + 1)] for j in range(2) for k in range(2)]
    # ll_rl ll_rr lr_rl lr_rr
    for j in range(2):
        for k in range(2):
            if d[j + 2 * k] == min(d):
                return [k == 1, j == 0]


def order_nn(Z, D):
    """Order leaves via knn algorithm."""
    for i in range(len(Z)):
        ordered = is_ordered(Z, D, i)
        for k in range(2):
            if not ordered[k]:
                swap_children(Z, i, k)
    return Z


dist_condensed = pdist(matrix, metric=distance)
cluster = order_nn(cluster, dist_condensed)

dend = dendrogram(cluster, no_plot=True)

# create list of genes without deleted genes (with zero expression in all samples)
calculated_genes = [gene for gene in all_genes if gene not in list_zero_genes]

calculated_genes_dict = {}
for i, gene in enumerate(calculated_genes):
    calculated_genes_dict[i] = {'gene': (gene)}

cluster = {
    'linkage': cluster.tolist(),
    'gene_symbols': calculated_genes_dict,
    'order': dend['leaves'],
    'zero_gene_symbols': list_zero_genes,
    'zero_sample_ids': [int(args.sampleids[index]) for index in sample_zero]
}

if args.output:
    with open(args.output, 'w') as f:
        json.dump(cluster, f)
else:
    print(json.dumps({'cluster': cluster}, separators=(',', ':')))

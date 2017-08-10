#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering of genes."""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, zscore
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from resolwe_runtime_utils import error, warning

from clustering_leaf_ordering import knn, optimal, simulated_annealing


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Hierarchical clustering of genes')
    parser.add_argument('-f', '--sample_files', nargs='+', help='Sample files', required=True)
    parser.add_argument('-s', '--sampleids', nargs='+', default=[], help='Sample IDs', type=int)
    parser.add_argument('-g', '--genes', nargs='+', default=[], help='Subset of gene labels')
    parser.add_argument('-d', '--distance', default='euclidean', help='Distance metric')
    parser.add_argument('-l', '--linkage', default='average', help='Linkage method')
    parser.add_argument('-o', '--ordering', default=None, help='Ordering method (knn, optimal, sa)')
    parser.add_argument('-k', '--n_keep', help='Keep _ states in optimal leaf ordering', type=int)
    parser.add_argument('-r', '--n_trials', help='Number of trials in simulated annealing',
                        type=int)
    parser.add_argument('-t', '--log2', action='store_true', help='Log2 transformation')
    parser.add_argument('-n', '--normalization', default=None, help='Normalization')
    parser.add_argument('--output', help='Output JSON filename')
    return parser.parse_args()


def get_expression(fname, sep='\t', gene_set=[]):
    """Read expressions from file and return only expressions of genes in gene_set."""
    df = pd.read_csv(fname, sep=sep, header=0, index_col=0, compression='gzip')
    df.index = df.index.map(str)
    if not gene_set:
        return df
    intersection = [gene for gene in gene_set if gene in df.index]
    return df.loc[intersection]


def get_expressions(fnames, sep='\t', gene_set=[]):
    """Read expressions from files.

    Return only expressions of genes that are listed in all samples and in the
    gene_set.

    """
    dfs = [get_expression(fname, sep=sep, gene_set=gene_set) for fname in fnames]
    df = pd.concat(dfs, axis=1)
    return df.dropna()


def transform(expressions, log2=False, const=1.0, normalization=None, ddof=1):
    """Compute log2 and normalize expression values.

    Parameters:
    - const: an additive constant used in computation of log2
    - normalization: None or 'z-score'
    - ddof: degrees of freedom used in computation of z-scores

    """
    if log2:
        try:
            expressions = expressions.applymap(lambda x: np.log2(x + const))
        except:
            msg = 'Cannot apply log2 to expression values.'
            print(error(msg))
            raise ValueError(msg)
    if normalization:
        if normalization == 'z-score':
            try:
                expressions = expressions.apply(lambda x: zscore(x, ddof=ddof), axis=1)
            except:
                msg = 'Cannot compute Z-scores.'
                print(error(msg))
                raise ValueError(msg)
        else:
            msg = 'Unknown normalization type {}.'.format(normalization)
            print(error(msg))
            raise ValueError(msg)
    return expressions


def get_distance_metric(distance_metric):
    """Get distance metric."""
    if distance_metric == 'spearman':
        return lambda x, y: 0.0 if x == y else 1.0 - spearmanr(x, y).correlation
    elif distance_metric in ['pearson', 'correlation']:
        return lambda x, y: 0.0 if x == y else 1.0 - pearsonr(x, y)[0]
    return distance_metric


def get_zero_genes(expressions):
    """Get genes which have zero expression in all samples."""
    return expressions[(expressions == 0.0).all(axis=1)].index.tolist()


def get_const_genes(expressions):
    """Get genes with constant expression profile across samples."""
    return expressions[expressions.apply(lambda x: min(x) == max(x), axis=1)].index.tolist()


def get_zero_samples(expressions):
    """Get indices of samples which have zero expression in all genes."""
    return [i for i, v in enumerate((expressions == 0.0).all()) if v]


def remove_zero_genes(expressions):
    """Remove genes with zero expression profile across samples."""
    return expressions.loc[(expressions != 0.0).all(axis=1)]


def remove_const_genes(expressions):
    """Remove genes with constant expression profile across samples."""
    return expressions[expressions.apply(lambda x: min(x) != max(x), axis=1)]


def get_clustering(expressions,
                   distance_metric='euclidean',
                   linkage_method='average',
                   ordering_method=None,
                   n_keep=None,
                   n_trials=1000):
    """Compute linkage, order, and produce a dendrogram."""
    if len(expressions.index) < 2:
        return np.array([]), {'leaves': list(expressions.index)}
    try:
        distances = pdist(np.array(expressions), metric=distance_metric)
        if np.isnan(distances).any():
            distances = np.nan_to_num(distances, copy=False)
            warning('Distances between some genes were undefined and were set to zero.')
    except:
        msg = 'Cannot compute distances between genes.'
        print(error(msg))
        raise ValueError(msg)
    try:
        link = linkage(y=distances, method=linkage_method)
    except:
        msg = 'Cannot compute linkage.'
        print(error(msg))
        raise ValueError(msg)
    if ordering_method:
        if ordering_method == 'knn':
            link = knn(link, distances)
        elif ordering_method == 'optimal':
            link = optimal(link, distances, n_keep)
        elif ordering_method == 'sa':
            link = simulated_annealing(link, distances, n_trials)
        else:
            msg = 'Unknown ordering method {}'.format(ordering_method)
            print(error(msg))
            raise ValueError(msg)
    try:
        dend = dendrogram(link, no_plot=True)
    except:
        msg = 'Cannot compute dendrogram.'
        print(error(msg))
        raise ValueError(msg)
    return link, dend


def output_json(result=dict(), fname=None):
    """Print json if fname=None else write json to file 'fname'."""
    if fname:
        with open(fname, 'w') as f:
            json.dump(result, f)
    else:
        print(json.dumps({'cluster': result}, separators=(',', ':')))


def main():
    """Compute gene hierarchical clustering."""
    args = parse_args()
    expressions = get_expressions(fnames=args.sample_files, gene_set=args.genes)
    expressions = transform(expressions, log2=args.log2, normalization=args.normalization)
    zero_genes = get_zero_genes(expressions)
    zero_samples = get_zero_samples(expressions)
    linkage, dendrogram = get_clustering(
        expressions,
        distance_metric=get_distance_metric(args.distance),
        linkage_method=args.linkage,
        ordering_method=args.ordering,
        n_keep=args.n_keep,
        n_trials=args.n_trials
    )
    result = {
        'linkage': linkage.tolist(),
        'gene_symbols': {i: {'gene': gene} for i, gene in enumerate(expressions.index)},
        'order': dendrogram['leaves'],
        'zero_gene_symbols': zero_genes,
        'missing_gene_symbols': list(set(args.genes).difference(set(expressions.index))),
        'zero_sample_ids': [args.sampleids[sample] for sample in zero_samples]
    }
    output_json(result, args.output)


# For a large number of genes the default recursion limit is hit during computation of dendrogram
sys.setrecursionlimit(10000)
main()

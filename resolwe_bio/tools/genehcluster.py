#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering of genes."""

import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, zscore
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

import resdk
from resolwe_runtime_utils import error, warning


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Hierarchical clustering of genes')
    parser.add_argument('-f', '--sample-files', nargs='+', help='Sample files', required=True)
    parser.add_argument('-n', '--sample-names', nargs='+', help='Sample names', required=True)
    parser.add_argument('-s', '--source', help='Source', required=True)
    parser.add_argument('-p', '--species', help='Species', required=True)
    parser.add_argument('-g', '--gene-labels', nargs='+', default=[], help='Subset of gene labels')
    parser.add_argument('-t', '--log2', action='store_true', help='Log2 transformation')
    parser.add_argument('-z', '--z-score', action='store_true', help='Z-score normalization')
    parser.add_argument('-r', '--remove-const', action='store_true', help='Remove genes with constant expression')
    parser.add_argument('-d', '--distance-metric', default='euclidean', help='Distance metric')
    parser.add_argument('-l', '--linkage-method', default='average', help='Linkage method')
    parser.add_argument('-o', '--order', action='store_true', help='Optimal ordering')
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

    Return only expressions of genes that are listed in all samples and in gene_set.

    """
    dfs = [get_expression(fname, sep=sep, gene_set=gene_set) for fname in fnames]
    inner = pd.concat(dfs, axis=1, join='inner')
    outer = pd.concat(dfs, axis=1, join='outer', sort=True)
    if gene_set:
        excluded = sorted(set(gene_set).difference(set(inner.index)))
    else:
        excluded = sorted(outer.index.difference(inner.index))
    return inner, excluded


def transform(expressions, log2=False, const=1.0, z_score=False, ddof=1):
    """Compute log2 and normalize expression values.

    Parameters:
    - log2: use log2(x+const) transformation
    - const: an additive constant used in computation of log2
    - z_score: use Z-score normalization
    - ddof: degrees of freedom used in computation of Z-score

    """
    if log2:
        expressions = expressions.applymap(lambda x: np.log2(x + const))
        if expressions.isnull().values.any():
            msg = 'Cannot apply log2 to expression values.'
            print(error(msg))
            raise ValueError(msg)
    if z_score:
        expressions = expressions.apply(lambda x: zscore(x, ddof=ddof), axis=1, result_type='broadcast')
        expressions.fillna(value=0.0, inplace=True)
    return expressions


def get_distance_metric(distance_metric):
    """Get distance metric."""
    if distance_metric == 'spearman':
        return lambda x, y: 1.0 - spearmanr(x, y).correlation
    elif distance_metric == 'pearson':
        return 'correlation'
    return distance_metric


def is_const(values):
    """Return True, if all values are approximately equal, otherwise return False."""
    mn = np.min(values)
    mx = np.max(values)
    if mn + mx == 0.0:
        return mn == mx
    else:
        return (mx - mn) / abs(mx + mn) < 1.0E-6


def remove_const_genes(expressions):
    """Remove genes with constant expression profile across samples."""
    matches = expressions.apply(lambda x: not is_const(x), axis=1)
    return expressions.loc[matches], matches.values.tolist()


def get_clustering(
        expressions,
        distance_metric='euclidean',
        linkage_method='average',
        order=False):
    """Compute linkage, order, and produce a dendrogram."""
    try:
        link = linkage(y=expressions, method=linkage_method, optimal_ordering=order)
    except:
        msg = 'Cannot compute linkage.'
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


def set_error(msg):
    """Print error message and raise ValueError."""
    print(error(msg))
    raise ValueError(msg)


def get_gene_names(feature_ids, source, species):
    """Map feature IDs to gene names."""
    res = resdk.Resolwe()
    features = res.feature.filter(feature_id=feature_ids, source=source, species=species)
    return [feature.name for feature in features]


def main():
    """Compute gene hierarchical clustering."""
    args = parse_args()

    if len(args.sample_files) != len(args.sample_names):
        msg = 'The number of sample files does not match the number of sample names.'
        set_error(msg)

    if len(args.gene_labels) == 1:
        msg = 'Select at least two genes to compute hierarchical clustering of genes.'
        set_error(msg)

    if len(args.sample_files) == 1 and args.distance_metric != 'euclidean':
        msg = ('Select at least two samples to compute hierarchical clustering of genes with '
               'correlation distance metric or use Euclidean distance metric.')
        set_error(msg)

    expressions, excluded = get_expressions(fnames=args.sample_files, gene_set=args.gene_labels)

    if len(expressions.index) == 0:
        if not args.gene_labels:
            msg = 'The selected samples do not have any common genes.'
        else:
            msg = 'None of the selected genes are present in all samples.'
        set_error(msg)

    if len(expressions.index) == 1 and args.distance_metric != 'euclidean':
        if not args.gene_labels:
            msg = ('The selected samples contain only one common gene ({}). At least two common '
                   'genes are required to compute hierarchical clustering of genes with '
                   'correlation distance metric. Select a different set of samples or use Euclidean '
                   'distance metric.'.format(get_gene_names(list(expressions.index), args.source, args.species)[0]))
        else:
            msg = ('Only one of the selected genes ({}) is present in all samples but at least two '
                   'such genes are required to compute hierarchical clustering of genes with '
                   'correlation distance metric. Select more genes or use Euclidean distance '
                   'metric.'.format(get_gene_names(list(expressions.index), args.source, args.species)[0]))
        set_error(msg)

    expressions = transform(expressions, log2=args.log2, z_score=args.z_score)

    if args.remove_const:
        expressions, matches = remove_const_genes(expressions)
        gene_names = get_gene_names(list(expressions.index), args.source, args.species)
        if len(expressions.index) == 0:
            msg = ('All of the selected genes have constant expression across samples. '
                   'Hierarchical clustering of genes cannot be computed.')
            set_error(msg)
        if len(expressions.index) == 1:
            msg = ('Only one of the selected genes ({}) has a non-constant expression across '
                   'samples. However, hierarchical clustering of genes cannot be computed with '
                   'just one gene.'.format(gene_names[0]))
            set_error(msg)
        removed = [name for i, name in enumerate(gene_names) if not matches[i]]
        suffix = '' if len(removed) <= 3 else ', ...'
        if removed:
            removed_names = get_gene_names(removed[:3], args.source, args.species)
            msg = ('{} of the selected genes ({}) have constant expression across samples. '
                   'Those genes are excluded from the computation of hierarchical clustering of '
                   'genes with correlation distance '
                   'metric.'.format(len(removed), ', '.join(removed_names) + suffix))
            print(warning(msg))
    else:
        matches = [True] * len(expressions.index)

    suffix = '' if len(excluded) <= 3 else ', ...'
    if excluded:
        excluded_names = get_gene_names(excluded[:3], args.source, args.species)
    if len(excluded) == 1:
        if not args.gene_labels:
            msg = ('Gene {} is present in some but not all of the selected samples. This '
                   'gene is excluded from the computation of hierarchical clustering of '
                   'genes.'.format(len(excluded), ', '.join(excluded_names)))
        else:
            msg = ('{} of the selected genes ({}) is missing in at least one of the selected '
                   'samples. This gene is excluded from the computation of hierarchical '
                   'clustering of genes.'.format(len(excluded), ', '.join(excluded_names)))
        print(warning(msg))
    if len(excluded) > 1:
        if not args.gene_labels:
            msg = ('{} genes ({}) are present in some but not all of the selected samples. Those '
                   'genes are excluded from the computation of hierarchical clustering of '
                   'genes.'.format(len(excluded), ', '.join(excluded_names)))
        else:
            msg = ('{} of the selected genes ({}) are missing in at least one of the selected '
                   'samples. Those genes are excluded from the computation of hierarchical '
                   'clustering of genes.'.format(len(excluded), ', '.join(excluded_names)))
        print(warning(msg))

    linkage, dendrogram = get_clustering(
        expressions,
        distance_metric=get_distance_metric(args.distance_metric),
        linkage_method=args.linkage_method,
        order=args.order
    )

    result = {
        'gene_symbols': {i: {'gene': gene} for i, gene in enumerate(expressions.index)},
        'linkage': linkage.tolist(),
        'order': dendrogram['leaves'],
    }
    output_json(result, args.output)


main()

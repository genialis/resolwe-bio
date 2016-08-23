#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,redefined-outer-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Principal components analysis."""
from __future__ import absolute_import, division, print_function

import argparse
import json
import sys

import numpy as np  # pylint: disable=import-error
from sklearn.decomposition import PCA  # pylint: disable=import-error

import utils


parser = argparse.ArgumentParser(description="PCA")
parser.add_argument('samples', help="All samples (comma separated)")
parser.add_argument('sample_ids', help="Sample IDs (comma separated")
parser.add_argument('--genes', nargs='+', help='filter genes')
parser.add_argument('--filter', help="Filter genes with low expression", action="store_true")
args = parser.parse_args()


def isfloat(value):
    """Check if value is float."""
    try:
        float(value)
        return True
    except ValueError:
        return False


def isgzipped(f):
    """Check if file f is gzipped."""
    with open(f, 'rb') as rpkm_file:
        magic = rpkm_file.read(2)

    return magic == '\037\213'


samples = args.samples.split(',')
sample_ids = args.sample_ids.split(',')

if len(samples) != len(sample_ids):
    print('{"rc":"1"}')
    exit(1)

exp = []
allgenes = set()

for fname in samples:
    myopen = utils.gzopen if isgzipped(fname) else open

    with myopen(fname) as f:
        exp.append({gene_exp[0]: float(gene_exp[1]) for gene_exp in
                    (l.split('\t') for l in f) if len(gene_exp) == 2 and isfloat(gene_exp[1])})

        allgenes.update(exp[-1].keys())

if args.genes:
    allgenes = allgenes.intersection(args.genes)

# Default expression value is 0.
allgenes_array = np.array(list(allgenes))
exp = np.array([[genemap.get(g, 0.) for g in allgenes_array] for genemap in exp])

if args.filter:
    exp = np.transpose(exp)
    f_exp = exp[np.sum(exp, axis=1) > exp.shape[1]]
    exp = np.transpose(f_exp)

if exp.shape[1] == 0:
    print(json.dumps({
        'proc.warning': 'Filtering removed all PCA attributes.',
        'pca': {
            'flot': {
                'data': [[0, 0] for i in range(exp.shape[0])],
                'xlabel': 'PC 1',
                'ylabel': 'PC 2',
                'samples': sample_ids}}}, separators=(',', ':')))

    sys.exit(0)

pca = PCA(n_components=0.99, whiten=True)
transformed_data = pca.fit_transform(exp)

coordinates = [[t[0], t[1]] if len(t) > 1 else [t[0], 0] for t in transformed_data]


def component_top_factors(component):
    """Return top 10 absolute factors."""
    # 10x faster, but not supported in current numpy:
    #   abs_component = np.abs(component)
    #   unordered_ixs = np.argpartition(abs_component, -10)[-10:]
    #   ixs = unordered_ixs[np.argsort(abs_component[unordered_ixs])[::-1]]
    ixs = np.argsort(np.abs(component))[:-11:-1]
    return zip(allgenes_array[ixs].tolist(), component[ixs].tolist())


data = {
    'pca': {
        'flot': {
            'data': coordinates,
            'xlabel': 'PC 1',
            'ylabel': 'PC 2',
            'samples': sample_ids
        }
    }
}
if not any(np.isnan(pca.explained_variance_ratio_)):
    data['pca']['all_explained_variance_ratios'] = pca.explained_variance_ratio_.tolist()
    data['pca']['all_components'] = [component_top_factors(component) for component in pca.components_]
    data['pca']['explained_variance_ratios'] = data['pca']['all_explained_variance_ratios'][:10]
    data['pca']['components'] = data['pca']['all_components'][:10]

print(json.dumps(data, separators=(',', ':'), allow_nan=False))

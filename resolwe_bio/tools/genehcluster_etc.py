#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
# XXX: Refactor to a comand line tool and remove pylint disable
"""Hierarchical clustering of expression time courses."""
from __future__ import absolute_import, division, print_function

import argparse
import json

import numpy as np
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform

# Import before Orange to avoid namespace conflicts
import utils

from Orange.clustering.hierarchical import HierarchicalClustering, AVERAGE, SINGLE, COMPLETE
from Orange.core import SymMatrix


parser = argparse.ArgumentParser(description='Hierarchical clustering of expression time courses.')

parser.add_argument('etc_files', nargs='+', help='etc files')
parser.add_argument('-e', '--expids', nargs='+', default=[], help='experiment ids')
parser.add_argument('-g', '--genes', nargs='+', default=[], help='subset of gene ids')
parser.add_argument('-d', '--dstfunc', default='pearson', help='distance function')
parser.add_argument('-l', '--linkage', default='average', help='clustering linkage function')
parser.add_argument('--output', help='Output JSON file')

args = parser.parse_args()

distance_map = {
    'spearman': lambda m: 0.5 - spearmanr(m, axis=1)[0] / 2.,
    'pearson': lambda m: 0.5 - np.corrcoef(m) / 2.,
    'euclidean': lambda m: squareform(pdist(m, 'euclidean')),
    'manhattan': lambda m: squareform(pdist(m, 'cityblock')),
}

linkage_map = {
    'average': AVERAGE,
    'single': SINGLE,
    'complete': COMPLETE
}

if args.dstfunc not in distance_map:
    raise ValueError("Invalid distance function {}".format(args.dstfunc))

if args.linkage not in linkage_map:
    raise ValueError("Invalid clustering linkage function {}".format(args.linkage))

if not args.expids or len(args.expids) != len(args.etc_files):
    raise ValueError("Number of experiment ids must match the number of files")

etcs = []
timepoints = set()

# read data
for i, fname in enumerate(args.etc_files):
    etcjson = json.load(utils.gzopen(fname))
    tps = etcjson['etc']['timePoints']
    expid = args.expids[i]

    if not all(tps[i] <= tps[i + 1] for i in range(len(tps) - 1)):
        raise ValueError("Timepoints should be ordered")

    etc = {'genes': {}, 'experiment': expid, 'timePoints': np.array(tps)}
    timepoints.update(tps)

    for gene in args.genes:
        if gene in etcjson['etc']['genes']:
            etc['genes'][gene] = np.array(etcjson['etc']['genes'][gene])

    etcs.append(etc)

timepoints = np.array(sorted(timepoints))
series, info = [], []

# interpolate missing timepoints
for etc in etcs:
    if not np.array_equal(timepoints, etc['timePoints']):
        for gene, points in etc['genes'].iteritems():
            series.append(np.interp(timepoints, etc['timePoints'], points))
            info.append((gene, etc['experiment']))
    else:
        for gene, points in etc['genes'].iteritems():
            series.append(points)
            info.append((gene, etc['experiment']))

matrix = distance_map[args.dstfunc](np.array(series))
matrix[np.isnan(matrix)] = np.nanmax(matrix)
matrix[matrix < 0] = 0.
matrix = SymMatrix([list(x) for x in matrix])

clustering = HierarchicalClustering()
clustering.linkage = linkage_map[args.linkage]
clustering.overwrite_matrix = 1
root = clustering(matrix)


def dendrogram(cluster):
    """Generate dendrogram structure."""
    res = {}
    q = [[cluster, res], ]

    while len(q) > 0:
        old, new = q.pop(0)

        if old.branches:
            new['left'] = {}
            new['right'] = {}
            new['height'] = old.height
            q.append([old.left, new['left']])
            q.append([old.right, new['right']])
        else:
            new['height'] = old.height
            new['gene'] = info[old[0]][0]
            new['experiment'] = info[old[0]][1]

    return res


clustering = {'tree': dendrogram(root)}

if args.output:
    with open(args.output, 'w') as f:
        json.dump(clustering, f)
else:
    print(json.dumps({'clustering': clustering}, separators=(',', ':')))

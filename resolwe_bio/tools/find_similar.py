#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Finding genes with simmilar expressions."""
from __future__ import absolute_import, division, print_function

import argparse
import json
import math

import numpy as np  # pylint: disable=import-error
from scipy.stats import spearmanr  # pylint: disable=import-error
from scipy.stats.stats import pearsonr  # pylint: disable=import-error

import utils


parser = argparse.ArgumentParser(description='Finding genes with simmilar expressions.')

parser.add_argument('-g', '--gene', help='gene ID')
parser.add_argument('-d', '--dstfunc', default='pearson', help='distance function')
parser.add_argument('etc_file', help='gene expression file')

args = parser.parse_args()


def euclidian(x, y):
    """Compute eucludean distance."""
    return np.linalg.norm(np.array(x) - np.array(y))


def pearson(x, y):
    """Compute Pearson's correlation."""
    return pearsonr(x, y)[0]


def spearman(x, y):
    """Compute Spearman's rank."""
    return spearmanr(x, y)[0]


# 2nd argument: True if higher value means better score
distance_map = {
    'euclidean': [euclidian, False],
    'pearson': [pearson, True],
    'spearman': [spearman, True]
}

search_gene = args.gene
file_handler = utils.gzopen(args.etc_file)
expressions = json.load(file_handler)
file_handler.close()
search_f, rev_sort = distance_map[args.dstfunc]

if args.dstfunc not in distance_map:
    raise ValueError("Invalid distance function {}".format(args.dstfunc))

search_gene_expression = expressions['etc']['genes'][search_gene]

similarity = [{'gene': gene, 'distance': search_f(expressions['etc']['genes'][gene], search_gene_expression)}
              for gene in expressions['etc']['genes'] if gene != search_gene]


similarity = [value for value in similarity if not math.isnan(value['distance'])]
similarity = sorted(similarity, key=lambda x: x['distance'], reverse=rev_sort)
similarity = {'search gene': search_gene, 'similar genes': similarity}
print(json.dumps({'simgenes': similarity}, separators=(',', ':')))

#!/usr/bin/env python2
import json
import math
import numpy as np
import argparse
import utils

from scipy.stats import spearmanr
from scipy.stats.stats import pearsonr

parser = argparse.ArgumentParser(description='Finding genes with simmilar expressions.')

parser.add_argument('-g', '--gene', help='gene ID')
parser.add_argument('-d', '--dstfunc', default='pearson', help='distance function')
parser.add_argument('etc_file', help='gene expression file')

args = parser.parse_args()


def euclidian(x, y):
    return np.linalg.norm(np.array(x) - np.array(y))


def pearson(x, y):
    return pearsonr(x, y)[0]


def spearman(x, y):
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


similarity = filter(lambda x: not math.isnan(x['distance']), similarity)
similarity.sort(reverse=rev_sort, key=lambda x: x['distance'])
similarity = {'search gene': search_gene, 'similar genes': similarity}
print json.dumps({'simgenes': similarity}, separators=(',', ':'))

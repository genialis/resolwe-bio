#!/usr/bin/env python2
import argparse
import json
import sys
import random


parser = argparse.ArgumentParser(description="Parse Diff Exp output files")
parser.add_argument('n', help="Number of genes")
args = parser.parse_args()

gene = ['gene_{}'.format(n_genes) for n_genes in xrange(int(args.n))]

de = {
    'status': ["OK" for n_genes in xrange(0, int(args.n))],
    'gene': gene,
    'gene_id': gene,
    'p_value': [random.uniform(0.000001, 1) for n_genes in xrange(int(args.n))],
    'value_1': [random.uniform(0, 10) for n_genes in xrange(int(args.n))],
    'value_2': [random.uniform(0, 10) for n_genes in xrange(int(args.n))],
    'test_stat': [random.uniform(-10, 10) for n_genes in xrange(int(args.n))],
    'sample_1': ['control' for n in xrange(int(args.n))],
    'sample_2': ['case' for n in xrange(int(args.n))],
    'locus': ['chr7:52823164-52845080' for n in xrange(int(args.n))],
    'q_value': [random.uniform(0,1) for n in xrange(int(args.n))],
    'significant': ['yes' for n in xrange(int(args.n))],
    'log2(fold_change)': [random.uniform(-10, 10) for n_genes in xrange(int(args.n))],
    'test_id': gene,
}

json_data = json.dumps({'de_data': de}, separators=(',', ':'), allow_nan=False)

size = sys.getsizeof(json_data) / 1024.**2
print json.dumps({'size': size}, separators=(',', ':'), allow_nan=False)

print json_data

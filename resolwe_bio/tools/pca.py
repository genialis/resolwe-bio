from collections import defaultdict
import gzip
import os
import argparse
import json

import numpy as np
from sklearn.decomposition import PCA

parser = argparse.ArgumentParser(description="PCA")
parser.add_argument('samples', help="All samples (comma separated)")
parser.add_argument('sample_ids', help="Sample IDs (comma separated")
args = parser.parse_args()


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def isgzipped(f):
    with open(f, 'rb') as rpkm_file:
        magic = rpkm_file.read(2)

    return magic == '\037\213'


samples = args.samples.split(',')
sample_ids = args.sample_ids.split(',')

if len(samples) != len(sample_ids):
    print '{"rc":"1"}'
    exit(1)

exp = []
allgenes = set()

for fname in samples:
    myopen = gzip.open if isgzipped(fname) else open

    with myopen(fname) as f:
        exp.append({gene_exp[0]: float(gene_exp[1]) for gene_exp in
                    (l.split('\t') for l in f) if len(gene_exp) == 2
                    and isfloat(gene_exp[1])})

        allgenes.update(exp[-1].keys())

# Default expression value is 0.
exp = np.array([[genemap.get(g, 0.) for g in allgenes] for genemap in exp])

pca = PCA(n_components=2)
transformed_data = pca.fit_transform(exp)

data = {'pca': {'flot': {'data': transformed_data.tolist(), 'xlabel': 'PC 1',
        'ylabel': 'PC 2', 'samples': sample_ids}}}

print json.dumps(data, separators=(',', ':'))

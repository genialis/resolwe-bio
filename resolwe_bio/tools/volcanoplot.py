import argparse
import json
import os

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Compute coordinates for volcano plot.')
parser.add_argument('de_results', help='Differential Expression results')

args = parser.parse_args()

if not os.path.isfile(args.de_results):
    print json.dumps({'proc.error': 'Missing DE results file.'}, separators=(',', ':'))
    exit(1)

def is_gzipped(f):
    with open(f, 'rb') as de_file:
        magic = de_file.read(2)

    return magic == '\037\213'

if args.de_results.endswith(('.xls', 'xlsx')):
    de = pd.io.excel.read_excel(args.de_results, sheetname=0)
elif is_gzipped(args.de_results):
    de = pd.read_csv(args.de_results, compression='gzip', sep="\t")
else:
    de = pd.read_csv(args.de_results, sep="\t")

de = de.dropna()
header = list(de)
ids = list(de.ix[:,0])

# get FC data
if 'log2FoldChange' in header:
    x = np.array(de['log2FoldChange'])
elif 'avg_FC' in header:
    x = np.array(de['avg_FC'])
elif 'logFC' in header:
    x = np.array(de['logFC'])

# get FDR/pval data
if 'padj' in header:
    y = -np.log10(np.array(de['padj']))
elif 'FDR.DE' in header:
    y = -np.log10(np.array(de['fdr.de']))
elif 'FDR' in header:
    y = -np.log10(np.array(de['FDR']))

try:
    data = {'volcano_plot': {'flot': {'data': zip(x, y)}, 'xlabel': 'log2',
        'ylabel': '-log10(FDR)', 'id': ids}}
    print json.dumps(data, separators=(',', ':'))
except NameError:
    print json.dumps({'proc.error': 'FC and/or FDR/pval data is missing.'}, separators=(',', ':'))
    exit(1)

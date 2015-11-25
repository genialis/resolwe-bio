#!/usr/bin/env python2
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
    xlabel = 'log2FoldChange'
elif 'paired_avg_FC' in header:
    x = np.array(de['paired_avg_FC'])
    xlabel = 'paired_avg_FC'
elif 'logFC' in header:
    x = np.array(de['logFC'])
    xlabel = 'logFC'
elif 'log2(fold_change)' in header:
    x = np.array(de['log2(fold_change)'])
    xlabel = 'log2(fold_change)'

# get FDR/pval data
if 'ebays.pval' in header:
    y = -np.log10(np.array(de['ebays.pval']))
    ylabel = '-log10(' + 'ebays.pval' + ')'
elif 'padj' in header:
    y = -np.log10(np.array(de['padj']))
    ylabel = '-log10(' + 'padj' + ')'
elif 'FDR.DE' in header:
    y = -np.log10(np.array(de['fdr.de']))
    ylabel = '-log10(' + 'fdr.de' + ')'
elif 'FDR' in header:
    y = -np.log10(np.array(de['FDR']))
    ylabel = '-log10(' + 'FDR' + ')'
elif 'q_value' in header:
    y = -np.log10(np.array(de['q_value']))
    ylabel = '-log10(' + 'q_value' + ')'

y[y == np.inf] = np.amax(y[np.isfinite(y)])

try:
    data = {'volcano_plot': {'flot': {'data': zip(x, y)}, 'xlabel': xlabel,
        'ylabel': ylabel, 'id': ids}}
    print json.dumps(data, separators=(',', ':'))
except NameError:
    print json.dumps({'proc.error': 'FC and/or FDR/pval data is missing.'}, separators=(',', ':'), allow_nan=False)
    exit(1)

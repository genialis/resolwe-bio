#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Compute coordinates for volcano plot."""
from __future__ import absolute_import, division, print_function

import argparse
import csv
import json
import os

import numpy as np  # pylint: disable=import-error

import utils


parser = argparse.ArgumentParser(description='Compute coordinates for volcano plot.')
parser.add_argument('bayseq_results', help='baySeq results')

args = parser.parse_args()

if not os.path.isfile(args.bayseq_results):
    exit(1)


def is_gzipped(f):
    """Check if file f is gzipped."""
    with open(f, 'rb') as rpkm_file:
        magic = rpkm_file.read(2)

    return magic == '\037\213'


base, ext = os.path.splitext(args.bayseq_results)
delimiter = ';' if ext == '.csv' else '\t'
ids = []
case_rpkum_median = []
control_rpkum_median = []
fdr_de = []

myopen = utils.gzopen if is_gzipped(args.bayseq_results) else open

with myopen(args.bayseq_results) as csvfile:
    reader = csv.reader(csvfile, delimiter=delimiter)
    header = reader.next()

    case_rpkum_median_index = header.index('Case_RPKUM_Median')
    control_rpkum_median_index = header.index('Control_RPKUM_Median')
    fdr_de_index = header.index('FDR.DE')

    for row in reader:
        case_rpkum_median.append(float(row[case_rpkum_median_index]))
        control_rpkum_median.append(float(row[control_rpkum_median_index]))
        fdr_de.append(float(row[fdr_de_index]))
        ids.append(row[0])

case_rpkum_median = np.array(case_rpkum_median)
control_rpkum_median = np.array(control_rpkum_median)

fdr_de = np.array(fdr_de)

x = np.log2((case_rpkum_median + 1) / (control_rpkum_median + 1))

# Zero would yield Infinity
fdr_de[fdr_de == 0] = np.min(fdr_de[fdr_de > 0])
y = -np.log10(fdr_de)

data = {'volcano_plot': {'flot': {'data': zip(x, y)}, 'xlabel': 'log2', 'ylabel': '-log10(FDR)', 'id': ids}}

print(json.dumps(data, separators=(',', ':')))

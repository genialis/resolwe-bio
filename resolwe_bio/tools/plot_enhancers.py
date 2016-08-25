#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Parse coordinates for the hockey-stick plot."""
from __future__ import absolute_import, division, print_function

import argparse
import json

import pandas as pd  # pylint: disable=import-error


parser = argparse.ArgumentParser(description='Parse coordinates for the hockey-stick plot.')
parser.add_argument('input_data', help='Input Data')
parser.add_argument('-c', action='store_true', help='Substract control sample data')

args = parser.parse_args()

xlabel = 'Enhancers ranked by signal'
ylabel = 'Total reads per million'
labels = []

data = pd.read_csv(args.input_data, header=5, sep='\t')

for line in open(args.input_data):
    if 'Using cutoff' in line:
        cutoff = line.strip().split(' ')[3]
        break

if args.c:
    x_axis = data.iloc[:, 8][::-1]
    y_axis = data.iloc[:, 6] - data.iloc[:, 7]
else:
    x_axis = data.iloc[:, 7][::-1]
    y_axis = data.iloc[:, 6]

n_sup_enh, rows = data[data.isSuper == 1].shape

chr_pos = data.CHROM.map(str) + ":" + data.START.map(str) + '-' + data.STOP.map(str)

if len(x_axis) != len(y_axis):
    print(json.dumps({'proc.error': 'Scatter plot error. len(x_axis) != len(y_axis)'}, separators=(',', ':')))
    exit(1)

if len(labels) > 0 and len(labels) != len(x_axis):
    print(json.dumps({'proc.error': 'Scatter plot error. len(labels) != len(x_axis)'}, separators=(',', ':')))
    exit(1)

data = {
    'scatter_plot': {
        'points': {
            'x_axis': list(x_axis),
            'y_axis': list(y_axis),
            'items': labels
        },
        'annotations': [
            {'type': 'line', 'x1': 0, 'y1': float(cutoff), 'x2': max(x_axis), 'y2': float(cutoff)},
            {'type': 'line', 'x1': len(x_axis) - n_sup_enh, 'y1': 0, 'x2': len(x_axis) - n_sup_enh, 'y2': max(y_axis)}
        ],
        'meta': {
            'x_label': xlabel,
            'y_label': ylabel,
            'chr_pos': list(chr_pos),
            'text': 'Cutoff: {}'.format(cutoff)
        }
    }
}

print(json.dumps(data, separators=(',', ':')))

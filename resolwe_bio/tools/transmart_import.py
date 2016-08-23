#!/usr/bin/env python
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Import gene expressions and annotations from tranSMART."""
from __future__ import absolute_import, division, print_function

import argparse
import csv
import json
import os
import re

from resolwe_runtime_utils import export  # pylint: disable=import-error

from . import transmart_utils


parser = argparse.ArgumentParser(description='Import gene expressions and the '
                                             'corresponding annotations from tranSMART.')
parser.add_argument('expressions', type=str, help='gene expressions file')
parser.add_argument('--ann', type=str, help='sample annotation file')
parser.add_argument('--ann_ids', type=str, help='annotation ids file')
parser.add_argument('--anntree', type=str, help='annotation attribute tree')
parser.add_argument('--progress', type=float, default=0., help='start progress')
args = parser.parse_args()

if not os.path.isfile(args.expressions):
    print('{{"proc.error": "Gene expressions file {} does not exist"}}'.format(args.expressions))

if args.ann and not os.path.isfile(args.ann):
    print('{{"proc.error": "Sample annotation file {} does not exist"}}'.format(args.ann))

if args.ann_ids and not os.path.isfile(args.ann_ids):
    print('{{"proc.error": "Annotation ids file {} does not exist"}}'.format(args.ann_ids))

if args.ann and not args.anntree:
    print('{{"proc.error": "Annotation attribute tree must be given with annotations"')

if args.anntree and not os.path.isfile(args.anntree):
    print('{{"proc.error": "Annotation attribute tree file {} does not exist"}}'.format(args.anntree))

var_samples, var_template = None, None

if args.ann:
    with open(args.ann, 'rb') as annfile:
        with open(args.anntree, 'rb') as treefile:
            with open(args.ann_ids, 'rb') as ann_ids_file:
                var_samples, var_template = transmart_utils.format_annotations(annfile, treefile, ann_ids_file)

with open(args.expressions, 'rb') as csvfile:
    exprs = csv.reader(csvfile, delimiter='\t', quotechar='"')
    header = exprs.next()
    exprs = zip(*list(exprs))

nsamples = len(exprs)
progress_step = nsamples / 10.
progress_milestone = progress_step
rartifact = re.compile('^X[0-9]')
gene_ids = [g[1:] if rartifact.match(g) else g for g in exprs[0]]  # pylint: disable=unsubscriptable-object

for i in range(1, nsamples):
    sample_id = header[i]
    if var_samples is not None and sample_id not in var_samples:
        continue

    fname = 'temp/{}.tab'.format(sample_id)

    with open(fname, 'wb') as tabfile:
        tabwriter = csv.writer(tabfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        tabwriter.writerow(['Gene', 'Expression'])
        tabwriter.writerows(zip(gene_ids, exprs[i]))  # pylint: disable=unsubscriptable-object

    d = {
        'process': 'upload-expression',
        'input': {
            'exp': os.path.basename(fname),
            'exp_type': 'Log2'
        }
    }

    if var_samples is not None:
        d['var_template'] = var_template
        d['var'] = var_samples[sample_id]

    print(export(fname))
    print('run {}'.format(json.dumps(d, separators=(',', ':'))))

    if i >= progress_milestone:
        print('{{"proc.progress": {}}}'.format(args.progress + (1. - args.progress) * i / nsamples))
        progress_milestone += progress_step

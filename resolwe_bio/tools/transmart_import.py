#!/usr/bin/env python
import argparse
import csv
import gzip
import json
import os
import re

import utils
import transmart_utils


parser = argparse.ArgumentParser(description='Import gene expressions and the '
                                             'corresponding annotations from tranSMART.')
parser.add_argument('expressions', type=str, help='gene expressions file')
parser.add_argument('--annotation', type=str, help='sample annotation file')
parser.add_argument('--progress', type=float, default=0., help='start progress')
args = parser.parse_args()

if not os.path.isfile(args.expressions):
    print '{{"proc.error": "Gene expressions file {} does not exist"}}'.format(args.expressions)

if args.annotation and not os.path.isfile(args.annotation):
    print '{{"proc.error": "Sample annotation file {} does not exist"}}'.format(args.annotation)

var_samples, var_template = None, None

if args.annotation:
    with open(args.annotation, 'rb') as csvfile:
        var_samples, var_template = transmart_utils.format_annotations(csvfile)

with open(args.expressions, 'rb') as csvfile:
    exprs = csv.reader(csvfile, delimiter='\t', quotechar='"')
    header = exprs.next()
    exprs = zip(*list(exprs))

nsamples = len(exprs)
progress_current = args.progress
progress_step = (1. - args.progress) / nsamples
rartifact = re.compile('^X[0-9]')
gene_ids = [g[1:] if rartifact.match(g) else g for g in exprs[0]]

for i in range(1, nsamples):
    sample_id = header[i]
    if var_samples is not None and sample_id not in var_samples:
        continue

    fname = 'temp/{}.tab'.format(sample_id)

    with open(fname, 'wb') as tabfile:
        tabwriter = csv.writer(tabfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        tabwriter.writerow(['Gene', 'Expression'])
        tabwriter.writerows(zip(gene_ids, exprs[i]))

    d = {
        'status': 'resolving',
        'processor_name': 'import:upload:expression',
        'input': {
            'exp': {
                'file': os.path.basename(fname),
                'file_temp': os.path.join(os.getcwd(), fname)
            },
            'exp_type': 'Log2'
        }
    }

    if var_samples is not None:
        d['var_template'] = var_template
        d['var'] = var_samples[sample_id]

    print 'run {}'.format(json.dumps(d, separators=(',', ':')))

    progress_current += progress_step
    print '{{"proc.progress": {}}}'.format(progress_current)

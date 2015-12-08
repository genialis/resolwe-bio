#!/usr/bin/env python
import argparse
import csv
import gzip
import json
import os
import xlrd

import utils


parser = argparse.ArgumentParser(description='Import gene expressions and the '
                                             'corresponding annotations from tranSMART.')
parser.add_argument('expressions', type=str, help='gene expressions file')
parser.add_argument('--annotation', type=str, help='sample annotation file')
parser.add_argument('--template', type=str, help='annotation template')
parser.add_argument('--progress', type=float, default=0., help='start progress')
args = parser.parse_args()

if not os.path.isfile(args.expressions):
    print '{{"proc.error": "Gene expressions file {} does not exist"}}'.format(args.expressions)

if args.annotation and not os.path.isfile(args.annotation):
    print '{{"proc.error": "Sample annotation file {} does not exist"}}'.format(args.annotation)

template = []
if args.template:
    try:
        template = json.loads(args.template.replace('&quot;', '"'))
    except ValueError:
        print '{"proc.error": "Error parsing annotation template"}'
        raise

sample_characteristics_ndx = -1

if args.annotation:
    wb_ann = xlrd.open_workbook(args.annotation)
    ws_ann = wb_ann.sheets()[0]
    attributes = ws_ann.row_values(0, 0)
    samples_ann = ws_ann.col_values(0, 1, ws_ann.nrows)

    try:
        sample_characteristics_ndx = attributes.index('Sample Characteristics Ch1')
    except ValueError:
        pass

has_annotation = sample_characteristics_ndx >= 0

wb = xlrd.open_workbook(args.expressions)
ws = wb.sheets()[0]

# Save Excel expressions as TAB
tabname, _ = os.path.splitext(args.expressions)
tabname += '.tab.gz'

with gzip.open(tabname, 'wb') as tabfile:
    tabwriter = csv.writer(tabfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    tabwriter.writerows(ws.row_values(row, 0, ws.ncols) for row in xrange(ws.nrows))

os.makedirs('temp')
print '{"expset": {"file": "%s", "refs": ["temp"]}}' % tabname

genes = ws.col_values(0, 1, ws.nrows)
samples = ws.row_values(0, 1)

sample_ndx = 1
progress_current = args.progress
progress_step = (1. - args.progress) / len(samples)

for sample in samples:
    exp = ws.col_values(sample_ndx, 1, ws.nrows)
    sample_ndx += 1

    fname = 'temp/{}.tab'.format(sample)

    with open(fname, 'wb') as tabfile:
        tabwriter = csv.writer(tabfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        tabwriter.writerow(['Gene', 'Expression'])
        tabwriter.writerows(zip(genes, exp))

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

    if template:
        d['var_template'] = template

    if has_annotation:
        var = {}
        try:
            anndx = samples_ann.index(sample) + 1
            for key_val in ws_ann.cell(anndx, sample_characteristics_ndx).value.split(' : '):
                vals = key_val.split(': ')

                if vals[0] in ('age', 'plateid'):
                    vals[1] = int(vals[1])

                var[utils.escape_mongokey(vals[0])] = vals[1]

        except ValueError:
            pass

        d['var'] = var

    print 'run {}'.format(json.dumps(d, separators=(',', ':')))

    progress_current += progress_step
    print '{{"proc.progress": {}}}'.format(progress_current)

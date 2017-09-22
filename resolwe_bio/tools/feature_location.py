#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Extract feature locations."""
from __future__ import absolute_import, division, print_function

import json
import argparse
import re

from six import iteritems

import utils


parser = argparse.ArgumentParser(description='Extract feature locations')
parser.add_argument('--annotation', help='Annotation file')
parser.add_argument('--feature_type', help='Feature type')
parser.add_argument('--id_type', help='Feature ID')
parser.add_argument('--summarize_exons', action='store_true', default=False)
args = parser.parse_args()

gene_id_re = re.compile(r'gene_id "([\w\-\.\$]*)"')
transcript_id_re = re.compile(r'transcript_id "([\w\-\.\$]*)"')
id_re = re.compile(r'ID=([\w\-\.\$]*)')
parent_re = re.compile(r'Parent=([\w\-\.\$]*)')


def _search(regex, string):
    match = regex.search(string)
    return match.group(1) if match else 'N/A'


def get_gene_id(ids):
    """Get gene id."""
    return _search(gene_id_re, ids)


def get_transcript_id(ids):
    """Get transcript id."""
    return _search(transcript_id_re, ids)


def get_id(ids):
    """Get id."""
    return _search(id_re, ids)


def get_parent_id(ids):
    """Get parent id."""
    return _search(id_re, ids)


feature_locations = {}

with open(args.annotation) as f:
    n_rows = sum(1 for _ in f)

count, nrows = 0, float(n_rows)
count_span = nrows / 100.
count_threshold = count_span

with open(args.annotation) as annotation:
    for l in annotation:
        a = l.split('\t')
        if len(a) < 8:
            continue
        if a[2] == args.feature_type:
            chromosome = a[0]
            start = a[3]
            end = a[4]

            if args.id_type == 'ID':
                feature_id = get_id(l)
            elif args.id_type == 'Parent':
                feature_id = get_parent_id(l)
            elif args.id_type == 'gene_id':
                feature_id = get_gene_id(l)
            elif args.id_type == 'transcript_id':
                feature_id = get_transcript_id(l)

            if args.summarize_exons:
                if feature_id not in feature_locations:
                    feature_locations[feature_id] = {'chr': chromosome, 'str': [], 'end': []}

                feature_locations[feature_id]['str'].append(start)
                feature_locations[feature_id]['end'].append(end)

            else:
                feature_locations[feature_id] = {'chr': chromosome, 'str': start, 'end': end}

        count += 1
        if count > count_threshold:
            print('{{"proc.progress": {}}}'.format(round(count / nrows, 2)))
            count_threshold += count_span

    if args.summarize_exons:
        for feature_id in feature_locations:
            feature_locations[feature_id]['str'] = str(min(map(int, feature_locations[feature_id]['str'])))
            feature_locations[feature_id]['end'] = str(max(map(int, feature_locations[feature_id]['end'])))

print(json.dumps({'feature_location': feature_locations}, separators=(',', ':')))

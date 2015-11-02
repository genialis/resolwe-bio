#!/usr/bin/env python2
import gzip
import json
import os
import sys


if len(sys.argv) != 2:
    print '{"rc":"1"}'
    exit(1)

fname = sys.argv[1]

if not os.path.isfile(fname):
    print '{"rc":"1"}'
    exit(1)


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def escape_geneid(key):
    return key.replace('$', u'\uff04').replace('.', u'\uff0e')


with gzip.open(fname, 'rb') as f:
    # Split lines by tabs
    # Ignore lines without a number in second column
    # Build a dictionary of gene-expression pairs
    exp = {'genes': {escape_geneid(gene_exp[0]): float(gene_exp[1]) for
                     gene_exp in (l.split('\t') for l in f) if
                     len(gene_exp) == 2 and isfloat(gene_exp[1])}}

print '{"exp_json":%s}' % json.dumps(exp, separators=(',', ':'))

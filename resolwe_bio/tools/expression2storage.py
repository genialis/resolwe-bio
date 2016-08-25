#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Save expressions to storage."""
from __future__ import absolute_import, division, print_function

import json
import os
import sys

import utils


if len(sys.argv) != 2:
    print('{"rc":"1"}')
    exit(1)

fname = sys.argv[1]

if not os.path.isfile(fname):
    print('{"rc":"1"}')
    exit(1)


def isfloat(value):
    """Check if value is float."""
    try:
        float(value)
        return True
    except ValueError:
        return False

with utils.gzopen(fname) as f:
    # Split lines by tabs
    # Ignore lines without a number in second column
    # Build a dictionary of gene-expression pairs
    exp = {'genes': {utils.escape_mongokey(gene_exp[0]): float(gene_exp[1]) for
                     gene_exp in (l.split('\t') for l in f) if
                     len(gene_exp) == 2 and isfloat(gene_exp[1])}}

print('{"exp_json":%s}' % json.dumps(exp, separators=(',', ':')))

#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Creates a GFF file from a Repeater Masker tabular file."""
from __future__ import absolute_import, division, print_function

import argparse
import csv
import sys


RMSK_IDX = {
    'bin': 0,
    'swScore': 1,
    'milliDiv': 2,
    'milliDel': 3,
    'milliIns': 4,
    'genoName': 5,
    'genoStart': 6,
    'genoEnd': 7,
    'genoLeft': 8,
    'strand': 9,
    'repName': 10,
    'repClass': 11,
    'repFamily': 12,
    'repStart': 13,
    'repEnd': 14,
    'repLeft': 15,
    'id': 16
}

parser = argparse.ArgumentParser(description='Creates a GFF file from a Repeater Masker tabular file.')
parser.add_argument('--rmsk', dest='rmsk_file', help='Repeater Masker tabular file')
args = parser.parse_args()

with open(args.rmsk_file, 'r') as f:
    rdr = csv.reader(f, delimiter='\t')

    print('##gff-version 3')

    for r in rdr:
        entry = (
            r[RMSK_IDX['genoName']],
            'RepeatMasker',
            r[RMSK_IDX['repClass']],
            int(r[RMSK_IDX['genoStart']]) + 1,   # RepeaterMasker's out file is 0-based, GFF3 is 1-based
            int(r[RMSK_IDX['genoEnd']]) + 1,
            r[RMSK_IDX['swScore']],
            r[RMSK_IDX['strand']],
            '.',
            'Name=%s' % r[RMSK_IDX['repName']]
        )
        print('\t'.join(map(str, entry)))

sys.exit(0)

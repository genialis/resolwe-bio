#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Rename cuffnorm output tables headers."""
from __future__ import absolute_import, division, print_function

import argparse
from operator import itemgetter

parser = argparse.ArgumentParser(description="Rename cuffnorm output tables headers")
parser.add_argument('replicates', help="Define replicates and groups")
parser.add_argument('sample_names', help="Sample names")
args = parser.parse_args()

replicates = args.replicates.split(",")
sample_names = args.sample_names.split(",")

samples = []
# This sort must produce the same order as cufflinks_sample_groups.py
for x in sorted(zip(replicates, sample_names), key=itemgetter(0)):
    samples.append('{}_Group{}'.format(x[1], x[0]))


def replace_labels(table_name, sample_list):
    """Replace original table headers which contain group names with sample names.

    :param str table_name: cuffnorm output table files
    :param sample_list: list of sample names
    :type sample_list: list

    """
    with open(table_name, 'r') as genes_table:
        genes_table.readline()
        file_lines = genes_table.readlines()

    with open(table_name, 'w') as genes_table:
        genes_table.write("tracking_id" + '\t')
        genes_table.write('\t'.join(sample_list) + '\n')
        genes_table.writelines(file_lines)


tables = [
    'cds.count_table',
    'cds.fpkm_table',
    'genes.count_table',
    'genes.fpkm_table',
    'isoforms.count_table',
    'isoforms.fpkm_table',
    'tss_groups.count_table',
    'tss_groups.fpkm_table',
]

for table in tables:
    replace_labels(table, samples)

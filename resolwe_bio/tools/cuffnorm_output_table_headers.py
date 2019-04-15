#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Rename cuffnorm output tables headers."""

import argparse
from operator import itemgetter

parser = argparse.ArgumentParser(description="Rename cuffnorm output tables headers")
parser.add_argument('replicates', help="Define replicates and groups")
parser.add_argument('sample_names', help="Sample names")
parser.add_argument('file_names', help="File names")
args = parser.parse_args()

replicates = args.replicates.split(",")
sample_names = args.sample_names.split(",")
file_names = args.file_names.split(",")

reps_samples = dict(zip(sample_names, replicates))
samples_files = dict(zip(file_names, sample_names))

files_ids = {}
with open('samples.table', 'r') as samples_table:
    samples_table.readline()
    lines = samples_table.readlines()
    for line in lines:
        values = line.split('\t')
        files_ids[values[0]] = values[1]


def replace_labels(table_name):
    """Replace original cuffnorm table headers with sample names.

    :param str table_name: cuffnorm output table files

    """
    with open(table_name, 'r') as genes_table:
        header = genes_table.readline().split('\t')
        file_lines = genes_table.readlines()

    # Replace cuffnorm automatic header values with sample names
    for i, label in enumerate(header[1:]):
        sample = samples_files[files_ids[label.strip('\n')]]
        header[i + 1] = '{}_Group{}'.format(sample, reps_samples[sample]).replace(" ", "_")

    with open(table_name, 'w') as genes_table:
        genes_table.write('\t'.join(header) + '\n')
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
    replace_labels(table)

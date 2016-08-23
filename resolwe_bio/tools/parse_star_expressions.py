#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Parse STAR expression (ReadsPerGene) file."""
import argparse
import csv


parser = argparse.ArgumentParser(description="Parse STAR expression (ReadsPerGene) file.")
parser.add_argument('exp_file', help="STAR output file (ReadsPerGene.tab).")
parser.add_argument('prefix', help="Output file prefix.")
args = parser.parse_args()

with open(args.exp_file, 'r') as exp_file:
    exp_reader = csv.reader(exp_file, delimiter="\t")

    rc_unstranded = open('{}_rc_unstranded.tab'.format(args.prefix), 'w')
    rc_stranded = open('{}_rc_stranded.tab'.format(args.prefix), 'w')
    rc_stranded_reverse = open('{}_rc_stranded_reverse.tab'.format(args.prefix), 'w')

    rc_unstranded = csv.writer(rc_unstranded, delimiter='\t')
    rc_stranded = csv.writer(rc_stranded, delimiter='\t')
    rc_stranded_reverse = csv.writer(rc_stranded_reverse, delimiter='\t')

    rc_unstranded.writerow(["Gene", "Expression"])
    rc_stranded.writerow(["Gene", "Expression"])
    rc_stranded_reverse.writerow(["Gene", "Expression"])

    for row in exp_reader:
        rc_unstranded.writerow([row[0], row[1]])
        rc_stranded.writerow([row[0], row[2]])
        rc_stranded_reverse.writerow([row[0], row[3]])

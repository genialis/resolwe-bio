#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Create gene expression profiles."""
from __future__ import absolute_import, division, print_function

import argparse
import logging
import multiprocessing
import os

import biox  # pylint: disable=import-error


parser = argparse.ArgumentParser(description='Create gene expression profiles.')
parser.add_argument('genome_file', help='genome fasta file')
parser.add_argument('gff_file', help='gff3 file')
parser.add_argument('length', type=int, help='read length')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')

args = parser.parse_args()


if not os.path.isfile(args.genome_file):
    raise ValueError("File {} does not exist".format(args.genome_file))

genome_name, genome_ext = os.path.splitext(args.genome_file)

if genome_ext != '.fasta':
    raise ValueError("Expected .fasta file got {}".format(args.genome_file))

if not os.path.isfile(args.gff_file):
    raise ValueError("File {} does not exist".format(args.gff_file))

if os.path.splitext(args.gff_file)[1] not in ['.gtf', '.gff3']:
    raise ValueError("Expected .gff file got {}".format(args.gff_file))

if args.verbose:
    biox.utils.verbosity(logging.INFO)


fasta_file = '{}_mappability_{}.fasta'.format(genome_name, args.length)
biox.data.prepare_fasta_mapability(args.genome_file, fasta_file, args.length)
print('{"proc.progress":0.3}')

b = biox.map.Bowtie()
b.set_m(1)
b.enable_v(0)
b.set_mode_fasta()
b.set_processors(min(multiprocessing.cpu_count(), 2))

b.map(
    '{}_index'.format(genome_name), fasta_file,
    '{}_mappability_{}'.format(genome_name, args.length), verbose=args.verbose)
print('{"proc.progress":0.7}')

f_gtf = open('foo.gtf', 'w')
parents = {}
with open(args.gff_file) as f:
    for line in f:
        if line.strip().startswith('#'):
            continue
        vals = line.split('\t')
        attrlist = vals[-1].strip().split(';')
        attrs = {}
        for att in attrlist:
            key, val = att.split('=')
            attrs[key] = val

        if vals[2] == 'CDS' or vals[2] == 'exon':
            gene_id = attrs['Parent']
            f_gtf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                vals[0], '', 'exon', vals[3], vals[4], '.', vals[6], '.',
                'gene_id "{}"; gene_name "{}"; gene_description "{}"; gene_type "{}";'.format(
                    gene_id, parents[gene_id]['name'],
                    parents[gene_id]['description'], parents[gene_id]['gene_type'])))
        else:
            parents[attrs['ID']] = {'gene_type': vals[2]}
            parents[attrs['ID']]['name'] = attrs['Name'] if 'Name' in attrs else attrs['ID']
            parents[attrs['ID']]['description'] = attrs['Note'] if 'Note' in attrs else ''

f_gtf.close()
print('{"proc.progress":0.75}')

results = biox.expression.gene_expression('foo.gtf', '{}_mappability_{}.bam'.format(genome_name, args.length))
gtf = biox.data.Gtf('foo.gtf')
coding_len = {}
gene_len = {}
for gene_id, gene in gtf.genes.iteritems():
    coding = 0
    for feature in gene.features:
        if feature.type != 'exon':
            continue
        coding += (feature.stop - feature.start + 1)
    coding_len[gene_id] = coding
    gene_len[gene_id] = gene.stop - gene.start + 1

print('{"proc.progress":0.85}')

f = open('{}_mappability_{}.tab'.format(genome_name, args.length), "wt")
header = ['gene_id', 'coverage', 'coding_len', 'gene_len', 'mapability']
f.write('\t'.join(header) + '\n')
gene_ids = results.keys()
gene_ids.sort()
for gene_id in gene_ids:
    coverage = results[gene_id]
    row = [gene_id, coverage, coding_len[gene_id], gene_len[gene_id], '%.3f' % (float(coverage) / coding_len[gene_id])]
    f.write('\t'.join(str(x) for x in row) + '\n')
f.close()

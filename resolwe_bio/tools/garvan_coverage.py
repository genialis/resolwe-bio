#!/usr/bin/env pypy
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Exon Coverage Report (Garvan Institute)."""
from __future__ import absolute_import, division, print_function

import argparse
import bisect
import gzip
import itertools
import json
import re
import sys

from collections import defaultdict
from operator import itemgetter
from subprocess import PIPE, Popen


parser = argparse.ArgumentParser(description='Exon Coverage Report (Garvan Institute)')
parser.add_argument('filter', type=int, help='coverage threshold filter')
parser.add_argument('bam', type=str, help='mapping file')
parser.add_argument('gtf', type=str, help='gene model file')
parser.add_argument('genes', nargs='+', help='list of genes')
parser.add_argument('--variants', type=str, help='variants file')

args = parser.parse_args()
genes = args.genes


def runcmd(command):
    """Run bash command."""
    p = Popen(args=command, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = p.communicate()

    if 'samtools: can\'t load index' in err:
        print('{"proc.error": "samtools: can\'t load index"}')
        exit(1)

    return out


def median(lst, lst_size=None):
    """Compute median."""
    lst = sorted(lst)
    if lst_size:
        n = lst_size
        n_diff = n - len(lst)
        i = (n - 1) // 2 - n_diff
        if i < 0:
            if i == -1 and not n % 2:
                return lst[0] / 2.
            return 0
    else:
        n = len(lst)
        i = (n - 1) // 2

    return lst[i] if n % 2 else (lst[i] + lst[i + 1]) / 2.


# read variants
variant_input = defaultdict(lambda: ([], []))
is_var_sorted = True
if args.variants:
    for l in open(args.variants):
        if l[0] == '#':
            continue
        v = l.split('\t', 3)
        variant_input[v[0]][0].append(int(v[1]))
        variant_input[v[0]][1].append(v[2])

        if len(variant_input[v[0]][0]) > 1 and variant_input[v[0]][0][-2] > variant_input[v[0]][0][-1]:
            is_var_sorted = False

if not is_var_sorted:
    for key, val in variant_input.items():
        val = zip(val[0], val[1])
        val = sorted(val)
        pos, names = zip(*val)
        variant_input[key] = (pos, names)

# read sequence names from BAM header
seqnames = set()
bam_header = runcmd('samtools view -H {}'.format(args.bam))
for line in bam_header.split('\n'):
    values = line.split()
    if len(values) > 1 and values[0] == '@SQ':
        seqnames.add(values[1].split(':')[1])

# read gene model
exons, missing_genes = [], set(genes)
parse_t = re.compile(r'transcript_id[ =]"([^"]*)"')
parse_g = re.compile(r'gene_name[ =]"([^"]*)"')
for l in open(args.gtf):
    gtf = l.split('\t')

    if len(gtf) == 9 and gtf[2] == 'exon' and gtf[0] in seqnames:
        g = parse_g.search(gtf[8])

        if g.group(1) not in genes:
            continue

        missing_genes.discard(g.group(1))

        t = parse_t.search(gtf[8])

        exons.append({
            'seq': gtf[0],
            'start': int(gtf[3]),
            'end': int(gtf[4]),
            'gene_name': g.group(1),
            'transcript_id': t.group(1)
        })

        if '{}-{}'.format(int(gtf[3]), int(gtf[4])) == '88928799-88929480':
            print(exons[-1])

print(json.dumps({'missing': {'genes': sorted(missing_genes)}}, separators=(',', ':')))

# compute coverage
try:
    trans_file = gzip.open('transcripts.txt.gz', 'wb')
    exons_file = gzip.open('exons.txt.gz', 'wb')
    variants_file = gzip.open('variants.txt.gz', 'wb')

    trans_file.write('\t'.join([
        'Sequence Name', 'Location', 'Gene', 'Transcript', 'Coverage Total',
        'Coverage Max', 'Coverage Mean', 'Coverage Median', 'Bases All',
        'Bases Covered (>{0}x)', 'Bases Covered (>{0}x) Ratio', 'Exons All',
        'Exons Covered (>{0}x)', 'Exons Covered (>{0}x) Ratio', 'Variants All',
        'Variants Covered (>{0}x)', 'Variants Covered (>{0}x) Ratio']).format(args.filter) + '\n')

    exons_file.write('\t'.join([
        'Sequence Name', 'Location', 'Gene', 'Exon', 'Transcript', 'Coverage Total',
        'Coverage Max', 'Coverage Min', 'Coverage Mean', 'Coverage Median',
        'Bases All', 'Bases Covered (>{0}x)', 'Bases Covered (>{0}x) Ratio',
        'Variants All', 'Variants Covered (>{0}x)', 'Variants Covered (>{0}x) Ratio']).format(args.filter) + '\n')

    variants_file.write('\t'.join([
        'Sequence Name', 'Location', 'Gene', 'Exon', 'Transcript', 'Coverage',
        'X-Ref', 'Above Filter']) + '\n')

    count, nexons = 0, float(len(exons))
    count_span = nexons / 100.
    count_threshold = count_span
    exon_count = 1
    exons.sort(key=itemgetter('gene_name', 'transcript_id'))
    for gene_name, gene_exons in itertools.groupby(exons, key=itemgetter('gene_name')):
        gene_exons = list(gene_exons)
        gene_start = min(gene_exons, key=itemgetter('start'))['start']
        gene_end = max(gene_exons, key=itemgetter('end'))['end']
        gene_seq = gene_exons[0]['seq']

        var_pos = variant_input.get(gene_seq, ([], []))[0]
        var_names = variant_input.get(gene_seq, ([], []))[1]

        # call samtools for gene coverage
        gene_coverage = runcmd('samtools depth -r {}:{}-{} {}'.format(gene_seq, gene_start, gene_end, args.bam))
        # if not gene_coverage:
        #     # hack to handle no gene coverage
        #     gene_coverage = '{}\t-1\t0'.format(gene_seq)

        pos, coverage = [], []
        for line in gene_coverage.splitlines():
            vals = line.split('\t')
            pos.append(int(vals[1]))
            coverage.append(int(vals[2]))

        for transcript_id, trans_exons in itertools.groupby(gene_exons, key=itemgetter('transcript_id')):
            trans_coverage_max, trans_coverage_total, trans_bases_all, trans_bases_above_filter = 0, 0, 0, 0
            exons_all, exons_above_filter = 0, 0
            trans_start, trans_end = sys.maxsize, 0
            trans_variants_all, trans_variants_above_filter = 0, 0
            trans_coverage = []

            for e in trans_exons:
                # samtools returns ordered coverage positions
                start_ndx = bisect.bisect_left(pos, e['start'])
                end_ndx = bisect.bisect(pos, e['end'], start_ndx)
                exon_coverage = coverage[start_ndx:end_ndx]
                exon_coverage_len = e['end'] - e['start'] + 1
                variants_above_filter = 0
                trans_coverage.extend(exon_coverage)

                if not exon_coverage:
                    # hack to handle no exon coverage
                    exon_coverage = [0]

                var_start_ndx = bisect.bisect_left(var_pos, e['start'])
                var_end_ndx = bisect.bisect(var_pos, e['end'], var_start_ndx)

                # find variant coverage
                for v in range(var_start_ndx, var_end_ndx):
                    var_ndx = bisect.bisect_left(pos, var_pos[v])
                    vc = coverage[var_ndx] if len(pos) > var_ndx and pos[var_ndx] == var_pos[v] else 0
                    if vc > args.filter:
                        variants_above_filter += 1

                    variants_file.write(
                        gene_seq + '\t' +
                        str(var_pos[v]) + '\t' +
                        gene_name + '\t' +
                        str(exon_count) + '\t' +
                        transcript_id + '\t' +
                        str(vc) + '\t' +
                        var_names[v] + '\t' +
                        str(vc > args.filter) + '\n'
                    )

                coverage_total = sum(exon_coverage)
                coverage_max = max(exon_coverage)
                bases_above_filter = len([c for c in exon_coverage if c > args.filter])
                variants_all = var_end_ndx - var_start_ndx

                exons_file.write(
                    e['seq'] + '\t' +
                    str(e['start']) + '-' + str(e['end']) + '\t' +
                    e['gene_name'] + '\t' +
                    str(exon_count) + '\t' +
                    e['transcript_id'] + '\t' +
                    str(coverage_total) + '\t' +
                    str(coverage_max) + '\t' +
                    str(min(exon_coverage)) + '\t' +
                    str(coverage_total / float(exon_coverage_len)) + '\t' +
                    str(median(exon_coverage, exon_coverage_len)) + '\t' +
                    str(exon_coverage_len) + '\t' +
                    str(bases_above_filter) + '\t' +
                    str(bases_above_filter / float(exon_coverage_len)) + '\t' +
                    str(variants_all) + '\t' +
                    str(variants_above_filter) + '\t' +
                    str((variants_above_filter / float(variants_all)) if variants_all else 0.) + '\n'
                )
                exon_count += 1
                print(exon_count, float(len(exons)))
                trans_coverage_total += coverage_total
                trans_coverage_max = max(trans_coverage_max, coverage_max)
                trans_bases_all += exon_coverage_len
                trans_bases_above_filter += bases_above_filter
                exons_all += 1
                exons_above_filter += 1 if bases_above_filter > 0 else 0
                trans_start = min(trans_start, e['start'])
                trans_end = max(trans_end, e['end'])
                trans_variants_all += variants_all
                trans_variants_above_filter += variants_above_filter

                count += 1
                if count > count_threshold:
                    print('{{"proc.progress": {}}}'.format(round(count / nexons, 2)))
                    count_threshold += count_span

            trans_file.write(
                gene_seq + '\t' +
                str(trans_start) + '-' + str(trans_end) + '\t' +
                gene_name + '\t' +
                transcript_id + '\t' +
                str(trans_coverage_total) + '\t' +
                str(trans_coverage_max) + '\t' +
                str((trans_coverage_total / float(trans_bases_all)) if trans_bases_all else 0.) + '\t' +
                str(median(trans_coverage, trans_bases_all)) + '\t' +
                str(trans_bases_all) + '\t' +
                str(trans_bases_above_filter) + '\t' +
                str((trans_bases_above_filter / float(trans_bases_all)) if trans_bases_all else 0.) + '\t' +
                str(exons_all) + '\t' +
                str(exons_above_filter) + '\t' +
                str((exons_above_filter / float(exons_all)) if exons_all else 0.) + '\t' +
                str(trans_variants_all) + '\t' +
                str(trans_variants_above_filter) + '\t' +
                str((trans_variants_above_filter / float(trans_variants_all)) if
                    trans_variants_above_filter else 0.) + '\n'
            )
finally:
    trans_file.close()
    exons_file.close()
    variants_file.close()

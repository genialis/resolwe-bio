#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Create BEDGRAPH coverage given GFF3 annotations."""
from __future__ import absolute_import, division, print_function

import argparse
import csv
import re
import sys

import utils


parser = argparse.ArgumentParser(
    description='Create BEDGRAPH coverage file for a tab file w.r.t. given GFF3 annotations.')
parser.add_argument('--tab', dest='tab_file', help='Tab file')
parser.add_argument('--tab-coverage-col', dest='tab_col_val', help='Tab column with coverage value')
parser.add_argument('--gff3', dest='gff3_file', help='GFF3 file')
args = parser.parse_args()

# Fetch gene ids and their expressions from tab file
with utils.gzopen(args.tab_file) as f:
    rdr = csv.reader(f, delimiter='\t')
    rdr.next()  # skip header
    tab_vals = {row[0]: float(row[int(args.tab_col_val)]) for row in rdr}

genes = {}
# Fetch gene regions and chromosomes they belong to
with open(args.gff3_file, 'r') as f:
    rdr = csv.reader(f, delimiter='\t')
    gene_id_regex = re.compile(r'ID=([A-Za-z0-9_]+);')
    for i, row in enumerate(rdr):
        # skip GFF3 headers
        if row[0][0:2] == '##':
            continue
        # skip if not mRNA
        if row[2] != 'mRNA' or row[2] != 'transcript':
            continue

        gene_id = gene_id_regex.search(row[8])
        if gene_id is None:
            print("No gene id found in line %d" % i)
            sys.exit(1)
        else:
            gene_id = gene_id.groups()[0]

        region = (int(row[3]), int(row[4]))
        if gene_id in genes:
            genes[gene_id]['regions'].append(region)
        else:
            genes[gene_id] = {
                'regions': [region],
                'chr': row[0]
            }

# Create a bedgraph (still with overlapping regions)
bedgraph = []
for gid in set(genes.keys()) & set(tab_vals.keys()):  # pylint: disable=consider-iterating-dictionary
    if tab_vals[gid] == 0.0:
        continue
    expr_level = tab_vals[gid]
    for region in sorted(genes[gid]['regions']):
        bed_data = (genes[gid]['chr'], region[0], region[1], expr_level, gid)
        bedgraph.append(bed_data)
bedgraph.sort(key=lambda t: (t[0], int(t[1]), int(t[2])))

# Flatten regions that overlap
last_reg_idx = len(bedgraph) - 1
unique_regions = []
i = 0
while i <= last_reg_idx:
    overlap_reg_chr = bedgraph[i][0]
    overlap_reg_start = bedgraph[i][1]
    overlap_reg_end = bedgraph[i][2]

    # working index set - indices of those regions that overlap in some continuous region
    ws_idxs = [i]
    j = i + 1
    while j <= last_reg_idx and bedgraph[j][1] <= overlap_reg_end and bedgraph[j][0] == overlap_reg_chr:
        ws_idxs.append(j)
        overlap_reg_end = max(overlap_reg_end, bedgraph[j][2])
        i = j
        j += 1

    def which_start(point):
        """Return indices of regions with $point as starting point."""
        return [_idx for _idx in ws_idxs if bedgraph[_idx][1] == point]

    def which_end(point):
        """Return indices of regions with $point as ending point."""
        return [_idx for _idx in ws_idxs if bedgraph[_idx][2] == point]

    if len(ws_idxs) > 1:
        starts = [bedgraph[_idx][1] for _idx in ws_idxs]
        ends = [bedgraph[_idx][2] for _idx in ws_idxs]
        points = list(set(starts + ends))

        # active index set - regions that are active in the current working index set for a given point p
        # (that is, p lies in all regions whose indices are in the active set)
        active_idxs = []

        def active_expr_avg():
            """Compute average expression."""
            return sum(bedgraph[_idx][3] for _idx in active_idxs) / len(active_idxs)

        # Crunch overlapping regions into unique, flat parts.
        prev_p = None
        for p in sorted(points):
            start_idx = which_start(p)
            end_idx = which_end(p)
            n_start = len(start_idx)
            n_end = len(end_idx)

            # current point is a starting point of some region(s) A and an ending point of some other region(s) B
            #
            #     mandatory ASCII art:
            #
            # A:  --------------------
            # B:                     -------------------------
            # C:                     --------------------------------
            #     ^                  ^
            #   prev_p               p
            #
            #     |-----------------||                                     <---- created unique & flat regions
            #            flat_1      flat_2
            #
            if n_start > 0 and n_end > 0:
                unique_regions.append((overlap_reg_chr, prev_p, p - 1, active_expr_avg()))      # flat_1
                active_idxs = list(set(active_idxs + start_idx))
                unique_regions.append((overlap_reg_chr, p, p, active_expr_avg()))               # flat_2
                active_idxs = [idx for idx in active_idxs if idx not in end_idx]
                prev_p = p + 1

            # current point is a starting point of some region(s)
            #
            #     mandatory ASCII art:
            #
            # A:  ---------------------------
            # B:                   -----------------------
            # C:                   ------------------------------
            #     ^                ^
            #   prev_p             p
            #
            #     |---------------|
            #          flat_1
            #
            elif n_start > 0:
                if prev_p is not None:
                    unique_regions.append((overlap_reg_chr, prev_p, p - 1, active_expr_avg()))
                active_idxs = list(set(active_idxs + start_idx))
                prev_p = p

            # current point is an ending point of some region(s)
            #
            #     mandatory ASCII art:
            #
            # A:  ---------------------
            # B:        ---------------
            # C:                --------------------------
            #                   ^     ^
            #                prev_p   p
            #
            #                   |-----|
            #                    flat_1
            #
            elif n_end > 0:
                unique_regions.append((overlap_reg_chr, prev_p, p, active_expr_avg()))
                active_idxs = [idx for idx in active_idxs if idx not in end_idx]
                prev_p = p + 1
    else:
        # No overlapping for this feature; can be safely added to unique regions.
        unique_regions.append((overlap_reg_chr, overlap_reg_start, overlap_reg_end, bedgraph[i][3]))
    i += 1

print('\n'.join('\t'.join(map(str, bg_line[0:4])) for bg_line in unique_regions))
sys.exit(0)

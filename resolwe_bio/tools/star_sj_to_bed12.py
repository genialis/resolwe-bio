#!/usr/bin/python3
"""Recalculate from STAR SJ.out.tab file to BED12 format."""

import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('sj_file', help="STAR SJ.out.tab output file")
args = parser.parse_args()

# STAR SJ.out.tab file consist of following columns:
# column 1: chromosome
# column 2: first base of the intron (1-based)
# column 3: last base of the intron (1-based)
# column 4: strand (0: undefined, 1: +, 2: -)
# column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
# column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
# column 7: number of uniquely mapping reads crossing the junction
# column 8: number of multi-mapping reads crossing the junction
# column 9: maximum spliced alignment overhang
sj_file = pd.read_csv(args.sj_file, delimiter='\t', header=None)

# BED12 consists of 12 columns:
header = [
    'chromosome',
    'sj_start',
    'sj_end',
    'sj_name',
    'score',
    'strand',
    'thick_start',
    'thick_end',
    'item_rgb',
    'block_counts',
    'block_sizes',
    'block_starts',
]

bed_file = pd.DataFrame(index=sj_file.index, columns=header)
# 1: chromosome = first column from STAR SJ.out.tab
bed_file.loc[:, 'chromosome'] = sj_file.iloc[:, 0].values

# 2: SJ start (0-based) =
# (first base of the intron (1-based) - maximum spliced alignment overhang) -1 (to recalculate to
# 0 based system)
bed_file.loc[:, 'sj_start'] = (sj_file.iloc[:, 1]) - (sj_file.iloc[:, 8]) - 1

# 3: SJ end (0-based) =
# (last base of the intron (1-based) + maximum spliced alignment overhang)
bed_file.loc[:, 'sj_end'] = (sj_file.iloc[:, 2]) + (sj_file.iloc[:, 8])

# 4: SJ name
rows_num_length = len(str(len(sj_file.index)))
bed_file.loc[:, 'sj_name'] = (
    (sj_file.index + 1).astype(str).map(lambda x: 'JUNC0' + x.zfill(rows_num_length))
)

# 5: score = number of uniquely and multi mapping reads crossing the junction
bed_file.loc[:, 'score'] = (sj_file.iloc[:, 6].values + sj_file.iloc[:, 7].values)

# 6: strand =  0: '.' (undefined) , 1: '+', 2: '-
conditions = [
    sj_file.iloc[:, 3] == 0,
    sj_file.iloc[:, 3] == 1,
    sj_file.iloc[:, 3] == 2
]
choices_strand = ['.', '+', '-']
bed_file.loc[:, 'strand'] = np.select(conditions, choices_strand)

# 7: thick start is the same as SJ start
bed_file.loc[:, 'thick_start'] = (sj_file.iloc[:, 1]) - (sj_file.iloc[:, 8]) - 1

# 8: thick end is the same as SJ end
bed_file.loc[:, 'thick_end'] = (sj_file.iloc[:, 2]) + (sj_file.iloc[:, 8])

# 9: item RGB = 255,0,0 (red color) for '-' strand, 0,0,255 (blue color) for '+' strand
# and 0,0,0 (black) for undefined
choices_rgb = ['0,0,0', '0,0,255', '255,0,0']
bed_file.loc[:, 'item_rgb'] = np.select(conditions, choices_rgb)

# 10: block counts = 2
bed_file.loc[:, 'block_counts'] = '2'

# 11: block sizes = maximum spliced alignment overhang, maximum spliced alignment overhang
bed_file.loc[:, 'block_sizes'] = (
    (sj_file.iloc[:, 8]).astype(str) + ',' + (sj_file.iloc[:, 8]).astype(str)
)

# 12: block starts (a comma-separated list of block starts, relative to SJ start)
# = 0, (SJ end - SJ start + maximum spliced alignment overhang +1 )
bed_file.loc[:, 'block_starts'] = (
    '0' +    # first block allways starts at SJ start
    ',' +
    ((sj_file.iloc[:, 2]) - (sj_file.iloc[:, 1]) + (sj_file.iloc[:, 8]) + 1).astype(str)
)

bed_file.to_csv('junctions_unsorted.bed', sep='\t', index=False, header=False)

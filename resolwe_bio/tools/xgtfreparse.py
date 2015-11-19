import argparse
import os
import sys

import biox


parser = argparse.ArgumentParser(description='Create gene expression profiles.')
parser.add_argument('fasta', help='input fasta file')
parser.add_argument('out', help='output fasta file')

args = parser.parse_args()

mask = {}
mask['chr2'] = [(3016085, 3768655)]
mask['chrBF'] = [(64985, 72996)]
mask['chrR'] = [(42801, 78150)]

genome = {}

f = biox.data.Fasta(args.fasta)

while f.read():
    chr = f.id.split('Chromosome: ')
    chr = chr[1].split(' ')[0]
    chr = 'chr%s' % chr
    chr_len = len(f.sequence)
    if mask.get(chr, -1) != -1:
        for (start, stop) in mask[chr]:
            f.sequence = f.sequence[0:start] + 'N' * (stop - start + 1) + f.sequence[stop + 1:]
    assert(chr_len == len(f.sequence))
    genome[chr] = f.sequence

chrs = sorted(genome.keys())

with open(args.out, 'wt') as f:
    for chr in chrs:
        f.write('>%s\n' % chr)
        i = 0
        while (i * 60) < len(genome[chr]):
            f.write('%s\n' % genome[chr][i * 60:(i + 1) * 60])
            i += 1

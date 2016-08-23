#!/usr/bin/env python
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Compatibility check for genome with alignment/annotaion."""
from __future__ import absolute_import, division, print_function

import argparse
from Bio import SeqIO  # pylint: disable=import-error


parser = argparse.ArgumentParser(description='Compatibility check for genome with alignment/annotaion')
parser.add_argument('-g', '--genome', help='Genome names in FASTA format')
parser.add_argument('-a', '--gtf', help='Annotation fiel in GFF or GTF file')
parser.add_argument('-b', '--mapping', help='BAM file header')
args = parser.parse_args()


def parse_bam_header(mapping):
    """Parse BAM header."""
    seqs = {}
    with open(mapping, 'r') as handle:
        for line in handle:
            values = line.split()
            if values[0] == '@SQ':
                seq = values[1].split(':')[1]
                length = values[2].split(':')[1]
                seqs[seq] = int(length)
    return seqs


def parse_genome(genome):
    """Parse genome."""
    seqs = {}
    with open(genome, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[record.id] = len(record)
    return seqs


def parse_annotations(annotation):
    """Parse annotations."""
    seqs = {}
    with open(annotation, 'r') as handle:
        for line in handle:
            if line[0] == "#":
                continue
            line = line.split()
            seqs[line[0]] = max([int(line[4]), seqs.get(line[0], 0)])
    return seqs


def compatibility_g_m(gen_1, gen_2, name_1, name_2):
    """Check compatibility."""
    print("Checking compatibility of {} with {}.".format(name_1, name_2))
    r_code = 0
    for seq in gen_1:
        if seq not in gen_2:
            print("FAIL\t{} sequence not found in {}.".format(seq, name_2))
            r_code = 2
        elif gen_1[seq] != gen_2[seq]:
            print("FAIL\t{} sequence lengths missmatch"
                  "(lengths are {} and {}).".format(seq, gen_1[seq], gen_2[seq]))
            r_code = 2
    print()
    return r_code


def compatibility_g_a(gen, anot):
    """Check compatibility with annotations."""
    print("Checking compatibility of genome with annotation file")
    r_code = 0
    for seq in gen:
        if seq not in anot:
            print("WARN\t{} sequence not found in annotaion file".format(seq))
            r_code = 1
    for seq in anot:
        if seq not in gen:
            print("FAIL\t{} sequence in annotation "
                  "but not in genome.".format(seq))
            r_code = 2
        elif anot[seq] > gen[seq]:
            print("FAIL\tannotation interval on {} sequence is out of "
                  "reference range.".format(seq))
            r_code = 2
    print()
    return r_code


if __name__ == "__main__":
    glob_r_code = 0
    genome_lens = parse_genome(args.genome)

    if args.mapping:
        mapping_lens = parse_bam_header(args.mapping)
        r_new = compatibility_g_m(genome_lens, mapping_lens, "genome", "BAM mapping")
        glob_r_code = max(r_new, glob_r_code)
        r_new = compatibility_g_m(mapping_lens, genome_lens, "BAM mapping", "genome")
        glob_r_code = max(r_new, glob_r_code)

    if args.gtf:
        annot_intv = parse_annotations(args.gtf)
        r_new = compatibility_g_a(genome_lens, annot_intv)
        glob_r_code = max(r_new, glob_r_code)

    if glob_r_code == 0:
        print("DONE.\tNO problems found.")
    if glob_r_code == 1:
        print("DONE.\tCheck warnings!")
    if glob_r_code == 2:
        print("DONE.\tCheck FAILS!")

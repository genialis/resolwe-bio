#!/usr/bin/env python3
"""
Remove overlapping indels from the LoFreq-derived VCF files.

Adapted after original script developed by Andreas Wilm (wilma@gis.a-star.edu.sg),
Genome Institute of Singapore. The original script was published under the MIT license.
"""

import argparse
import gzip
import sys

from collections import namedtuple


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Remove overlaping indels.")
    parser.add_argument('in_vcf', help="Input VCF file.")
    parser.add_argument('out_vcf', help="Output (filterd) VCF file.")
    return parser.parse_args()


def write_var(var, out_file):
    """Write variant to output file."""
    var = var._replace(pos=str(var.pos))
    out_file.write('{}\n'.format('\t'.join(var)))


def vcf_line_to_var(line):
    """Parse vcf lines to variants."""
    VCFEntry = namedtuple('VCFEntry', ['chrom', 'pos', 'dbsnpid', 'ref', 'alt', 'qual', 'filter', 'info'])
    fields = line.rstrip().split('\t')[:8]
    e = VCFEntry._make(fields)
    return e._replace(pos=int(e.pos))


def af_from_var(var):
    """Parse AF value."""
    for f in var.info.split(';'):
        if f.startswith('AF='):
            return float(f[3:])
    return None


def qual_from_var(var):
    """Take care of missing values, int conversion and ties in comparisons."""
    if var.qual == ".":
        return sys.maxint
    else:
        # add AF to deal with ties
        return int(var.qual) + af_from_var(var)


def overlap(v1, v2):
    """Determine whether affected positions of two variants overlap."""
    v1_b = v1.pos + max(len(v1.ref), len(v1.alt))
    v2_b = v2.pos + max(len(v2.ref), len(v2.alt))
    return min(v1_b, v2_b) - max(v1.pos, v2.pos) > 0


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    with open(args.in_vcf) as in_vcf:
        with open(args.out_vcf, 'w') as out_vcf:
            pick_best_func = qual_from_var

            prev_vars = []
            for line in in_vcf:
                line = line.rstrip()

                # Write header lines to the output VCF file
                if line.startswith('#'):
                    out_vcf.write('{}\n'.format(line))
                    continue

                cur_var = vcf_line_to_var(line)
                if len(prev_vars):
                    if cur_var.chrom != prev_vars[-1].chrom or not overlap(prev_vars[-1], cur_var):
                        # pick highest qual/af from stack and empty stack
                        picked_var = sorted(prev_vars, key=lambda e: pick_best_func(e), reverse=True)[0]
                        write_var(picked_var, out_vcf)
                        prev_vars = []
                prev_vars.append(cur_var)

            # don't forget remaining ones
            picked_var = sorted(prev_vars, key=lambda e: pick_best_func(e), reverse=True)[0]
            write_var(picked_var, out_vcf)


if __name__ == "__main__":
    main()

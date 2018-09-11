#!/usr/bin/python3
"""Change chromosome names in BigWig header from ENSEMBL or NCBI to UCSC."""

import argparse
import os
import sys

import pyBigWig
from resolwe_runtime_utils import warning


MAPPINGS_DIR = '/opt/chrom_mappings/assets/'
MAPPINGS_FILES = {
    'Homo sapiens': [
        'GRCh38.p12_ensembl2UCSC.txt',
        'GRCh38.p12_NCBI2UCSC.txt',
    ],
    'Mus musculus': [
        'GRCm38.p6_ensembl2UCSC.txt',
        'GRCm38.p6_NCBI2UCSC.txt',
    ],
    'Rattus norvegicus': [
        'rn6_ensembl2UCSC.txt',
    ],
}


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--infile', help="Input bigWig file")
    parser.add_argument('--outfile', help="Output bigWig file name")
    parser.add_argument('--species', help="Species")
    return parser.parse_args()


def parse_mapping_file(filename):
    """Parse mapping file."""
    mappings = {}
    with open(filename) as handle:
        for line in handle:
            line = line.strip().split()
            if len(line) < 2 or line[1] == "":
                # Skip line, if it does not contain names of both sources.
                continue
            mappings[line[0]] = line[1]
    return mappings


def parse_mappings(species, infile, outfile):
    """Parse file with chromosome mappings."""
    mappings = dict()
    # if species doesn't have prepared mapping file the script should exit with status 0 and return BigWig file
    # with output name and warining
    if species not in MAPPINGS_FILES:
        msg = 'Chromosome mappings for Species "{}" are not supported.'.format(species)
        print(warning(msg))
        os.rename(infile, outfile)
        sys.exit(0)

    for basename in MAPPINGS_FILES[species]:
        filename = os.path.join(MAPPINGS_DIR, basename)
        mappings.update(parse_mapping_file(filename))
    return mappings


def create_new_header(infile, mappings, outfile):
    """Create new header in BigWig, with UCSC chromosome names."""
    with pyBigWig.open(infile) as bw:
        if set(bw.chroms().keys()).issubset(mappings.values()):
            # If chromosome names are already UCSC, just rename input file to output name.
            # Exit with status 0 since this is normal behavior.
            os.rename(infile, outfile)
            sys.exit(0)

        hdr = [(mappings[chrom], length) for chrom, length in bw.chroms().items() if chrom in mappings]

        if not hdr:
            msg = "Neither of the chromosomes in the input file has a valid UCSC pair. No mapping will be done."
            print(warning(msg))
            os.rename(infile, outfile)
            sys.exit(0)

        seq_num = 0
        with pyBigWig.open(outfile, 'w') as bw_output:
            bw_output.addHeader(hdr)
            for chrom, length in bw.chroms().items():
                ints = bw.intervals(chrom, 0, length)
                if ints and chrom in mappings:
                    bw_output.addEntries([mappings[chrom]] * len(ints),
                                         [x[0] for x in ints],
                                         ends=[x[1] for x in ints],
                                         values=[x[2] for x in ints])
                elif chrom not in mappings:
                    seq_num += 1
                    print('UCSC chromosome/conting mapping for {} is missing'.format(chrom))

        if seq_num > 0:
            print(warning("UCSC chromosome/conting mapping for {} sequence(s) is missing. "
                          "This sequence(s) will not be included in the bigWig file.".format(seq_num)))


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()
    mappings = parse_mappings(args.species.strip("'").capitalize(), args.infile, args.outfile)
    create_new_header(args.infile, mappings, args.outfile)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate input FASTQ file(s)."""

import argparse
import collections

import dnaio

from os.path import basename, isfile
from resolwe_runtime_utils import error
from dnaio.exceptions import FastqFormatError, FileFormatError


SUPPORTED_EXTENSIONS = (
    '.fastq',
    '.fastq.gz',
    '.fq',
    '.fq.gz',
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Validate single/paired-end FASTQ files.")
    parser.add_argument('-fq', nargs='+', type=str, required=True, help="FASTQ file(s).")
    parser.add_argument('-fq2', nargs='+', type=str, required=False, help="Second-in-pair FASTQ file(s).")
    return parser.parse_args()


def set_error(msg):
    """Print error message and raise ValueError."""
    print(error(msg))
    raise ValueError(msg)


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    if args.fq2:
        input_fastq = args.fq + args.fq2
    else:
        input_fastq = args.fq

    for fq_file in input_fastq:
        if not isfile(fq_file):
            set_error('Input file {} does not exist'.format(basename(fq_file)))
        if not fq_file.lower().endswith(SUPPORTED_EXTENSIONS):
            set_error('Unrecognized file name extension in file {}. '
                      'Supported file name extensions are {}.'.format(fq_file, SUPPORTED_EXTENSIONS))

    # Reduce the probability of uploading the FASTQ files with the same
    # content multiple times (as multiple lanes or mates).
    if len(set(input_fastq)) != len(input_fastq):
        seen_files = [item for item, count in collections.Counter(input_fastq).items() if count > 1]
        set_error('Non-unique input file names detected: {}.'.format(seen_files))

    if args.fq2 and len(args.fq) != len(args.fq2):
        set_error('The number of mate-pair files in split-lane samples must match. '
                  '{} and {} input files were given for the -fq and -fq2 inputs, '
                  'respectively.'.format(len(args.fq), len(args.fq2)))

    if args.fq2:
        for mate1, mate2 in zip(args.fq, args.fq2):
            try:
                paired_reads = dnaio.open(mate1, file2=mate2, fileformat='fastq')
                if not any(paired_reads):
                    set_error(
                        'Mate-pair files {} and {} contain no read sequences.'.format(basename(mate1), basename(mate2))
                    )
                else:
                    for read in paired_reads:
                        continue
                    print('Successfully validated mate-pair files {} and {}.'.format(basename(mate1), basename(mate2)))

            except (FastqFormatError, FileFormatError) as dnaio_error:
                set_error(
                    'Format error in mate-pairs {} and {}. {}'.format(
                        basename(mate1), basename(mate2), str(dnaio_error))
                )
    else:
        for fq in args.fq:
            try:
                reads = dnaio.open(fq, fileformat='fastq')
                if not any(reads):
                    set_error('Input file {} contains no read sequences.'.format(basename(fq)))
                else:
                    for read in reads:
                        continue
                    print('Successfully validated reads file {}.'.format(basename(fq)))

            except (FastqFormatError, FileFormatError) as dnaio_error:
                set_error('Error in file {}. {}'.format(basename(fq), str(dnaio_error)))


if __name__ == "__main__":
    main()

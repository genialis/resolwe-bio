#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""gtf2gff.py -- A script to convert GTF to GFF files.

A better replacement for gtf2gff.pl

Version 1.0
"""
from __future__ import with_statement
import argparse

__author__ = "Uli Koehler & Anton Smirnov"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__ = "Apache v2.0"


class GTFException(Exception):
    """GTFException."""

    pass


def gtf2gff(infileName, outfileName, program):
    """Transform gtf to gff."""
    with open(infileName, "r") as infile, open(outfileName, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            words = line.split("\t")
            if len(words) != 9:
                raise GTFException("Found %d columns instead of the expected 9 in line: '%s'" % (len(words), line))
            else:
                outfile.write("\t".join(words) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--program', help='The name of the program which generated the GTF file.', required=True)
    parser.add_argument('infile', help="The GTF input file.",)
    parser.add_argument('outfile', help="The GFF output file.", nargs="?")
    args = parser.parse_args()
    gtf2gff(args.infile, args.outfile, args.program)

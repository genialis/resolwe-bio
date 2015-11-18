#!/usr/bin/env python2
"""
   awk 'NR % 4 == 0' your.fastq | python %prog [options]

guess the encoding of a stream of qual lines.
"""
import sys
import optparse


def get_qual_range(qual_str):
    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)


def get_encodings_in_range(lowestChar, rmax):
    if lowestChar < 33:
        return "Unknown_to_low"
    elif lowestChar < 64:
        return "Sanger"
    elif lowestChar == 65:
        return "Illumina_old"
    elif lowestChar <= 126:
        return "Illumina_old"
    else:
        return "Unknown_to_high"


def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-n", dest="n", help="number of qual lines to test default:-1"
                 " means test until end of file",
                 type='int', default=-1)

    opts, args = p.parse_args()
    gmin, gmax = 200, 0
    valid = []
    for i, line in enumerate(sys.stdin):
        lmin, lmax = get_qual_range(line.rstrip())
        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
            valid = get_encodings_in_range(gmin, gmax)
        if opts.n > 0 and i > opts.n:
            print valid
            sys.exit()
        if valid == "Unknown_to_low":
            print "Sanger"
            sys.exit()
        if valid == "Sanger":
            print "Sanger"
            sys.exit()
    print valid


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()

import os
import argparse
import logging

import pysam

import biox


parser = argparse.ArgumentParser(description='Create bigwig files for JBrowse.')
parser.add_argument('bam_file', help='bam file')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
args = parser.parse_args()

name, ext = os.path.splitext(args.bam_file)

if ext != '.bam':
    raise ValueError()

if args.verbose:
    biox.utils.verbosity(logging.INFO)


def bam2wig(bam_filename, bw_filename, strand=None, position='span', scale=None):
    strand_dic = {1: '+', -1: '-', '1': '+', '-1': '-', '+': '+', '-': '-'}
    chrs_filename = bam_filename + ".chrs"
    bed_filename = bw_filename + ".bed"
    write_bam_chr(bam_filename, chrs_filename)
    strand_parameter = "-strand %s" % strand_dic[strand] if strand is not None else ''
    scale_parameter = "-scale %.5f" % float(scale) if scale is not None else ''  # values are multiplied (*) by scale

    command = "genomeCoverageBed -bg -ibam {} -g {} {} {} > {}".format(
        bam_filename, chrs_filename, strand_parameter, scale_parameter, bed_filename)
    output, error = biox.utils.cmd(command)
    if args.verbose:
        print output

    command = "bedGraphToBigWig {} {} {}".format(bed_filename, chrs_filename, bw_filename)
    output, error = biox.utils.cmd(command)
    if args.verbose:
        print output

    os.remove(bed_filename)
    os.remove(chrs_filename)


bw_file = args.bam_file.replace(".bam", ".bwn")
bw_file_forward = args.bam_file.replace(".bam", "_f.bwn")
bw_file_reverse = args.bam_file.replace(".bam", "_r.bwn")
if args.verbose:
    print "computing mapped reads..."

mapped_reads = reduce(
    lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:]))
                         for l in pysam.idxstats(args.bam_file)])
scale = 1 / float(mapped_reads) * 1e6  # multiply by million

if args.verbose:
    print args.bam_file, mapped_reads, scale

biox.expression.bam2wig(args.bam_file, bw_file, scale=scale)
biox.expression.bam2wig(args.bam_file, bw_file_forward, strand="+", scale=scale)
biox.expression.bam2wig(args.bam_file, bw_file_reverse, strand="-", scale=scale)

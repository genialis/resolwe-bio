import argparse
import inspect
import logging
import multiprocessing
import os

import biox

parser = argparse.ArgumentParser(description='BCM enhanced short read aligner.')

parser.add_argument('data_dir', help='data object folder')
parser.add_argument('index_file', help='reference genome index')
parser.add_argument('read_file', help='FastQ reads file to map')
parser.add_argument('result_file', help='output file name')

parser.add_argument('-u', dest='num_u', type=int, default=0, help='map all reads')
parser.add_argument('--mode', choices=['n', 'v'], default='n', help='(n) use qualities or (v) use mismatches')
parser.add_argument('-v', dest='num_v', type=int, default=2, help='allowed mismatches')
parser.add_argument('-n', dest='num_n', type=int, default=2, help='allowed mismatches')
parser.add_argument('-l', dest='num_l', type=int, default=28, help='seed length')

parser.add_argument('-k', dest='num_k', type=int, default=1, help='alignements per read to report')
parser.add_argument('--strata', action='store_true', help='report only one best stratum')

parser.add_argument('--trim3', action='store_true', help='trim unmapped reads from 3\'')
parser.add_argument('--trim3-nucl', type=int, default=2, help='nucleotides to trim')
parser.add_argument('--trim3-iter', type=int, default=5, help='trim iterations')

parser.add_argument('--verbose', action='store_true', help='verbose output')

args = parser.parse_args()

if args.verbose:
    biox.utils.verbosity(logging.INFO)

print "\n--- BCM aligner: PARAMETERS ---\n"
for d in dir(args):
    arg = getattr(args, d)
    if not inspect.ismethod(arg) and not d.startswith('_'):
        print "%s: %s" % (d, arg)

# decode bowtie parameters and start mapping
b = biox.map.Bowtie()

b.set_m(args.num_k)
# use 1 core if only 1 available
b.set_processors(min(multiprocessing.cpu_count(), 2))
b.set_quality("phred64-quals")

if args.num_u:
    b.enable_u(args.num_u)

if args.mode == 'v':
    b.enable_v(args.num_v)
else:
    b.enable_n(args.num_n)
    b.set_l(args.num_l)

if args.trim3:
    b.enable_trim3(args.trim3_iter, args.trim3_nucl)

b.enable_strata(args.strata)

# if mapping.filetype=="FASTA":
#     b.set_mode_fasta()

print "\n--- BCM aligner: ALIGNING ---\n"
b.map(args.index_file, args.read_file, os.path.join(args.data_dir, args.result_file), verbose=args.verbose)

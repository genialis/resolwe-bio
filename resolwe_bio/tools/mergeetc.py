#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Save ETC data to output file."""
import json
import argparse
import os
import itertools

from collections import defaultdict

import utils


parser = argparse.ArgumentParser(description='Save ETC data to output file.')
parser.add_argument('files', nargs='*', help='ETC files')
parser.add_argument('--experiments', nargs='+', help='experiment ids')
parser.add_argument('--genes', nargs='+', help='filter genes')
parser.add_argument('--out', help='output file')

args = parser.parse_args()

out_file = open(args.out, "w")

header = ["Gene"]
geneset = []
exp = defaultdict(list)

experiments = iter(args.experiments)

for etc in args.files:
    if not os.path.isfile(etc):
        exit(1)

    with utils.gzopen(etc) as f:
        etc_data = json.load(f)
        x = next(experiments)
        header = header + [x + ' - ' + tp + 'h' for tp in map(str, etc_data["etc"]["timePoints"])]
        gn = set([g for g in etc_data["etc"]["genes"]])

        geneset.append(gn)

        for g in gn:
            exp[g].append(etc_data["etc"]["genes"][g])

genes = set.intersection(*geneset)

if args.genes:
    genes = genes.intersection(args.genes)

genes = sorted(genes)

out_file.write("\t".join(header) + '\n')

if len(genes) == 0:
    out_file.write("No common genes in selected ETC experiments.")
    print('{"proc.warning":"No common genes in selected ETC experiments."}')

for g in genes:
    exps = list(itertools.chain(*exp[g]))
    out_file.write("{}\t{}\n".format(g, "\t".join(map(str, exps))))

out_file.close()

#!/usr/bin/env python2
"""Parse GAF file to gene sets."""
import argparse
import csv
import json

from collections import defaultdict

from resolwe_runtime_utils import export  # pylint: disable=import-error
from utils import gzopen  # pylint: disable=import-error


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Create gene sets from gaf file.")
    parser.add_argument('gaf_file', help="GAF file.")
    parser.add_argument('source', help="GeneID source.")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    terms = defaultdict(list)

    with gzopen(args.gaf_file) as gaf:
        for line in gaf:
            if line.startswith('!'):
                continue
            go_data = line.strip().split('\t')
            terms[go_data[4]].append(go_data[1])

    for go_term in terms:
        outfile = '{}.tab'.format(str(go_term).replace(':', ''))
        with open(outfile, "wt") as f:
            writer = csv.writer(f, delimiter=str('\t'), lineterminator='\n')
            for gene in terms[go_term]:
                writer.writerow([gene])

            print(export(outfile))
            d = {
                'process': 'upload-geneset',
                'input': {
                    'src': outfile,
                    'source': args.source
                }
            }

        print('run {}'.format(json.dumps(d, separators=(',', ':'))))

    print('{{"num_genesets":{}}}'.format(len(terms)))


if __name__ == "__main__":
    main()

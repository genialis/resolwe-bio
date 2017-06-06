#!/usr/bin/env python3

"""RNA-Seq for biologists."""

from argparse import ArgumentParser
from json import dumps
from sys import exit
from resdk import Resolwe

GENOMES_AND_ANNOTATIONS = {
    'hs': ('genome-hg19', 'annotation-hg19'),
    'mm': ('genome-mm10', 'annotation-mm10')
}
ADAPTERS = {
    'yes': 'adapters-illumina',
    'no': None
}


def parse_arguments():
    """Parse command line arguments."""
    parser = ArgumentParser(description='Run RNA-Seq for biologists')
    parser.add_argument('type', help='Reads data type')
    parser.add_argument('reads', help='Reads')
    parser.add_argument('genome_and_annotation', help='Genome and annotation')
    parser.add_argument('adapters', help='Adapters')
    parser.add_argument('minlen', type=int, help='Min length')
    parser.add_argument('trailing', type=int, help='Trailing quality')
    parser.add_argument('stranded', help='Stranded')
    return parser.parse_args()


def get_data_id(res, slug, type=None):
    """Throw an error if data object is missing."""
    data_objs = res.data.filter(slug=slug, type=type)
    if not data_objs:
        print("Server configuration error: no slug '{}' with type '{}'.".format(slug, type))
        exit()
    return data_objs[0].id


def main():
    """Run RNA-Seq workflow."""
    args = parse_arguments()
    data = {
        'input': {
            'reads': args.reads,
            'minlen': args.minlen,
            'trailing': args.trailing,
            'stranded': args.stranded
        }
    }
    if args.type.startswith('data:reads:fastq:single:'):
        data['process'] = 'workflow-rnaseq-single'
    if args.type.startswith('data:reads:fastq:paired:'):
        data['process'] = 'workflow-rnaseq-paired'

    gna = GENOMES_AND_ANNOTATIONS[args.genome_and_annotation]
    ada = ADAPTERS[args.adapters]
    res = Resolwe()
    data['input']['genome'] = get_data_id(res, gna[0], 'data:genome:fasta:')
    data['input']['annotation'] = get_data_id(res, gna[1], 'data:annotation:gtf:')
    data['input']['adapters'] = get_data_id(res, ada, 'data:seq:nucleotide:')
    print('run {}'.format(dumps(data, separators=(',', ':'))))


if __name__ == '__main__':
    main()

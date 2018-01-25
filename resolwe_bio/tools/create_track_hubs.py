#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name
"""Create UCSC track hub files."""

import argparse
import os


def get_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Create files for a UCSC track hub")
    parser.add_argument('--name', nargs='+', help="Sample names")
    parser.add_argument('--track_hub_name', help="Track hub name")
    parser.add_argument('--build', help="Genome build")
    parser.add_argument('--browser_url', help="Local UCSC genome browser URL")
    parser.add_argument('--email', help="User email")
    return parser.parse_args()


def main():
    """Create track hub files."""
    args = get_args()
    build = args.build
    name = args.track_hub_name
    email = args.email
    if not email:
        email = 'support@genialis.com'

    trackdb_file = os.path.join(build, 'trackDb.txt')
    genome_file = '.'.join([name, 'genomes', 'txt'])
    with open(genome_file, 'w') as outfile:
        outfile.write('genome {}\ntrackDB {}'.format(build, trackdb_file))

    hub_file = '{}.hub.txt'.format(name)
    with open(hub_file, 'w') as outfile:
        for label in ['hub', 'shortLabel', 'longLabel']:
            outfile.write(' '.join([label, name]) + '\n')
        sequence = '\n'.join([' '.join(['genomeFile', genome_file]), ' '.join(['email', email])])
        outfile.write(sequence)

    URL = 'http://taco-wiki.grid.bcm.edu/genialis/ucsc'
    with open(trackdb_file, 'w') as outfile:
        for sample in args.name:
            lines = [
                ' '.join(['track', sample]),
                ' '.join(['bigDataUrl', '/'.join([URL, name, build, '{}.bw'.format(sample)])]),
                ' '.join(['shortLabel', sample]),
                ' '.join(['longLabel', sample]),
                'type bigWig',
                'color 0,0,0'
            ]
            outfile.write('\n'.join(lines))
            outfile.write('\n\n')

    trackhub_url = '/'.join([URL, name, hub_file])

    print('/'.join([args.browser_url.strip('/'), 'cgi-bin/hgHubConnect?hubUrl=', trackhub_url]))


if __name__ == '__main__':
    main()

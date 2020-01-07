#!/usr/bin/env python3
"""Format sample metadata table for MultiQC."""
import argparse
import json


sample_summary_json = {
    'id': 'sample_info',
    'section_name': 'Sample Info',
    'plot_type': 'table',
    'file_format': 'json',
    'data': {}
}


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', nargs='+', help='Samples.')
    parser.add_argument('--species', nargs='+', help='Species.')
    parser.add_argument('--build', nargs='+', help='Build.')
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    for sample_name, species, build in zip(args.samples, args.species, args.build):
        if build not in ['rRNA', 'globin']:
            sample_summary_json['data'][sample_name] = {
                'Species': species,
                'Genome Build': build,
            }

    with open('sample_data_mqc.json', 'w') as out_file:
        json.dump(sample_summary_json, out_file)


if __name__ == "__main__":
    main()

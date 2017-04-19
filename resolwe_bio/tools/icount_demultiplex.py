#!/usr/bin/env python3
"""Run iCount demultiplexing."""
from subprocess import Popen, PIPE, DEVNULL

import argparse
import gzip
import json
import os

from resolwe_runtime_utils import export  # pylint: disable=import-error


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run iCount demultiplexing.")
    parser.add_argument('annotation_file', help="iCount sample annotation file.")
    parser.add_argument('reads_file', help="Multiplexed reads file (fastq.gz)")
    return parser.parse_args()


def _tsv_to_dict(table_file):
    """Parse *.tsv file into dict."""
    data = {}
    with gzip.open(table_file, 'rt') as tfile:
        header = next(tfile).strip().split('\t')
        column_picker = {label: index for index, label in enumerate(header)}
        for line in tfile:
            line = line.split('\t')
            barcode = line[column_picker["5' barcode"]]
            brcs_list = barcode.split(',')
            brc = '{}{}{}'.format(brcs_list[0], brcs_list[1].split('_')[0], brcs_list[2])
            data[brc] = {key: line[column_picker[key]] for key in column_picker}
    return data


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    demux_data = _tsv_to_dict(args.annotation_file)
    barcodes = list(demux_data.keys())

    adapters_temp = list(set([demux_data[barcode]["3' adapter"] for barcode in barcodes]))
    if len(adapters_temp) != 1:
        print('{"proc.error":"Only one adapter type is allowed per demultiplexing job."}')
        exit(1)

    adapter = ''.join([adpr.split('_')[0] for adpr in adapters_temp[0].split(',')])

    process = Popen(['iCount', 'demultiplex', args.reads_file, adapter] + barcodes,
                    stdout=PIPE,
                    stderr=DEVNULL)
    out, err = process.communicate()

    demux_files = [f for f in os.listdir('.') if f.startswith('demux')]

    for demux_file in demux_files:
        print(export(demux_file))
        d = {
            'process': 'upload-fastq-single',
            'input': {
                'src': [demux_file]
            }
        }

        print('run {}'.format(json.dumps(d, separators=(',', ':'))))


if __name__ == "__main__":
    main()

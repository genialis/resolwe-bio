#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name,redefined-outer-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""NGS reads demultiplexer."""

import argparse
import gzip
import json
import os
import subprocess
import sys

from six import iteritems

from resolwe_runtime_utils import export  # pylint: disable=import-error


parser = argparse.ArgumentParser(description='NGS reads demultiplexer.')

parser.add_argument('barcodes', help='barcodes file')
parser.add_argument('s', metavar='READS', nargs='?', help='file containing unpaired reads')
parser.add_argument('-1', metavar='READS-1', help='file containing upstream mates')
parser.add_argument('-2', metavar='READS-2', help='file containing downstream mates')
parser.add_argument('-m', '--mapping', help='barcode mapping file')
parser.add_argument('--progress-start', type=float, default=0., help="initial progress")
args = parser.parse_args()

if (not (args.s or (args.__dict__['1'] and args.__dict__['2']))
        or (args.s and args.__dict__['1'] and args.__dict__['2'])):
    sys.stderr.write('Give either unpaired reads or both paired read mates.')
    print()
    exit(1)

if args.s:
    reads1 = args.s
    reads2 = ''
else:
    reads1 = args.__dict__['1']
    reads2 = args.__dict__['2']
    if not os.path.isfile(reads2):
        sys.stderr.write('Reads file {} not found.'.format(reads2))
        print()
        exit(1)

if not os.path.isfile(reads1):
    sys.stderr.write('Reads file {} not found.'.format(reads1))
    print()
    exit(1)

if not os.path.isfile(args.barcodes):
    sys.stderr.write('Barcodes file {} not found.'.format(args.barcodes))
    print()
    exit(1)

if args.mapping and not os.path.isfile(args.mapping):
    sys.stderr.write('Barcode mapping file {} not found.'.format(args.mapping))
    print()
    exit(1)

pool_maps = {}


def isnum(number):
    """Check if number."""
    try:
        int(number)
        return True
    except ValueError:
        return False


barcode_length = 0

if args.mapping:
    with open(args.mapping, encoding='utf-8') as fd:
        for line in fd:
            line = line.rstrip()
            if not line:
                continue

            t = line.split('\t')
            barcode, filename = '', ''

            if len(t) == 2:
                barcode, filename = t[0:2]

            if len(t) > 2 and isnum(t[0]):
                barcode, filename = t[1:3]

            barcode, filename = barcode.strip(), filename.strip()

            if barcode and filename:
                pool_maps[barcode] = filename

                if barcode_length > 0 and barcode_length != len(barcode):
                    print('{"proc.error":"Barcodes should be of the same length."}')
                    exit(1)
                else:
                    barcode_length = len(barcode)

for bar, _map in iteritems(pool_maps):
    print('{}: {}'.format(bar, _map))


def read_multiplexed(reads1_file, reads2_file, barcodes_file, pool_maps, progress_start):
    """Parse multiplexed file."""
    pool_name = reads1_file.split('.')[0]

    def nicename(a):
        return a.replace('#', '').replace('  ', ' ').replace('/', ' ').replace(' ', '_')

    files, f1, f2, fbar = {}, None, None, None
    try:
        barcodes = set(pool_maps.keys())
        print("BARCODES: {}".format(barcodes))

        for barcode in barcodes:
            name = nicename(pool_maps[barcode])
            if reads2_file:
                filename = '{}_{}_{}_mate1.fq.gz'.format(pool_name, name, barcode)
                files[barcode] = gzip.open(filename, 'wb')

                filename = '{}_{}_{}_mate2.fq.gz'.format(pool_name, name, barcode)
                files[barcode + '2'] = gzip.open(filename, 'wb')

            else:
                filename = '{}_{}_{}.fq.gz'.format(pool_name, name, barcode)
                files[barcode] = gzip.open(filename, 'wb')

        if reads2_file:
            files['notmatched'] = gzip.open('Not_Matched_{}_mate1.fq.gz'.format(pool_name), 'wb')
            files['badquality'] = gzip.open('Bad_Quality_{}_mate1.fq.gz'.format(pool_name), 'wb')
            files['notmatched2'] = gzip.open('Not_Matched_{}_mate2.fq.gz'.format(pool_name), 'wb')
            files['badquality2'] = gzip.open('Bad_Quality_{}_mate2.fq.gz'.format(pool_name), 'wb')
        else:
            files['notmatched'] = gzip.open('Not_Matched_{}.fq.gz'.format(pool_name), 'wb')
            files['badquality'] = gzip.open('Bad_Quality_{}.fq.gz'.format(pool_name), 'wb')

        filenames = list(sorted(set(f.name for f in files.values())))

        p = subprocess.Popen(
            'gzip -dc {} | wc -l'.format(barcodes_file),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        numlines, err = p.communicate()

        if err:
            raise Exception(err)

        numlines = int(numlines)
        readid, matched, notmatched, badquality, skipped = 0, 0, 0, 0, 0
        print(json.dumps({'proc.progress': progress_start}))
        progress = progress_start
        progress_step = (0.9 - progress) / 20.
        progress_span = numlines / 20

        def save_results(matched, notmatched, badquality, skipped, total, progress):
            total = float(total)

            print(json.dumps({
                'matched': '{:,} reads ({:.2f} %)'.format(matched, 100 * matched / total),
                'notmatched': '{:,} reads ({:.2f} %)'.format(notmatched, 100 * notmatched / total),
                'badquality': '{:,} reads ({:.2f} %)'.format(badquality, 100 * badquality / total),
                'skipped': '{:,} reads ({:.2f} %)'.format(skipped, 100 * skipped / total),
                'proc.progress': progress
            }, separators=(',', ':')))

        f1 = gzip.GzipFile(reads1_file, 'r')
        fbar = gzip.GzipFile(barcodes_file, 'r')

        if reads2_file:
            f2 = gzip.GzipFile(reads2_file, 'r')

        while True:
            readid += 1
            r1 = f1.readline()
            if not r1:
                break
            r1 = r1.decode('utf-8').rstrip('\r').rstrip('\n').split('\t')
            if len(r1) != 11:
                print("SKIPPED: error in {} line in r1".format(readid))
                continue
            s1 = r1[-3].replace('.', 'N')
            p1 = r1[-1]

            rbar = fbar.readline()
            if not rbar:
                break
            rbar = rbar.decode('utf-8').rstrip('\r').rstrip('\n').split('\t')
            if len(rbar) != 11:
                print("SKIPPED: error in {} line in rbar".format(readid))
                continue
            sbar = rbar[-3].replace('.', 'N')[:barcode_length]
            pbar = rbar[-1]

            if reads2_file:
                r2 = f2.readline()
                if not r2:
                    break
                r2 = r2.decode('utf-8').rstrip('\r').rstrip('\n').split('\t')
                if len(r2) != 11:
                    print("SKIPPED: error in {} line in r2".format(readid))
                    continue
                s2 = r2[-3].replace('.', 'N')
                p2 = r2[-1]
            else:
                r2 = r1
                p2 = p1

            if r1[:7] == r2[:7] == rbar[:7] and p1 == p2 == pbar:
                idline = '@' + ':'.join(r1[:7]) + ' ' + sbar
                if p1 == '1' and p2 == '1':
                    if sbar in barcodes:
                        files[sbar].write(
                            (idline + '\n' + s1 + '\n' + '+' + '\n' + r1[-2] + '\n').encode('utf-8'))
                        if reads2_file:
                            files[sbar + '2'].write(
                                (idline + '\n' + s2 + '\n' + '+' + '\n' + r2[-2] + '\n').encode(
                                    'utf-8'))
                        matched += 1
                    else:
                        files['notmatched'].write(
                            (idline + '\n' + s1 + '\n' + '+' + '\n' + r1[-2] + '\n').encode('utf-8'))
                        if reads2_file:
                            files['notmatched2'].write(
                                (idline + '\n' + s2 + '\n' + '+' + '\n' + r2[-2] + '\n').encode(
                                    'utf-8'))
                        notmatched += 1
                else:
                    files['badquality'].write(
                        (idline + '\n' + s1 + '\n' + '+' + '\n' + r1[-2] + '\n').encode('utf-8'))
                    if reads2_file:
                        files['badquality2'].write(
                            (idline + '\n' + s2 + '\n' + '+' + '\n' + r2[-2] + '\n').encode('utf-8'))
                    badquality += 1
            else:
                print("SKIPPED: {}, p1: {}, p2: {}, pbar: {}".format(readid, p1, p2, pbar))
                print("{} ? {} ? {}".format(r1[:7], r2[:7], rbar[:7]))
                skipped += 1

            if readid % progress_span == 0:
                progress += progress_step
                save_results(matched, notmatched, badquality, skipped, readid, progress)

        save_results(matched, notmatched, badquality, skipped, readid, 0.9)

    finally:
        if f1:
            f1.close()
        if f2:
            f2.close()
        if fbar:
            fbar.close()

        for f in files:
            files[f].close()

    return filenames


filenames = read_multiplexed(reads1, reads2, args.barcodes, pool_maps, args.progress_start)

for name in filenames:
    if reads2:
        if name.endswith('_mate2.fq.gz'):
            continue

        name2 = name.replace('_mate1', '_mate2')
        print(export(name))
        print(export(name2))
        d = {
            'process': 'upload-fastq-paired',
            'input': {
                'src1': [name],
                'src2': [name2]
            }
        }
    else:
        print(export(name))
        d = {
            'process': 'upload-fastq-single',
            'input': {
                'src': [name]
            }
        }

    print('run {}'.format(json.dumps(d, separators=(',', ':'))))

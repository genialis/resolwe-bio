import argparse
import glob
import gzip
import json
import os
import shutil
import subprocess
import sys

from collections import defaultdict


parser = argparse.ArgumentParser(description='NGS reads demultiplexer.')

parser.add_argument('barcodes', help='barcodes file')
parser.add_argument('s', metavar='READS', nargs='?', help='file containing unpaired reads')
parser.add_argument('-1', metavar='READS-1', help='file containing upstream mates')
parser.add_argument('-2', metavar='READS-2', help='file containing downstream mates')
parser.add_argument('-m', '--mapping', help='barcode mapping file')
parser.add_argument('--progress-start', type=float, default=0., help="initial progress")
args = parser.parse_args()

if (not (args.s or (args.__dict__['1'] and args.__dict__['2'])) or
        (args.s and args.__dict__['1'] and args.__dict__['2'])):
    sys.stderr.write('Give either unpaired reads or both paired read mates.')
    print
    exit(1)

if args.s:
    reads1 = args.s
    reads2 = ''
else:
    reads1 = args.__dict__['1']
    reads2 = args.__dict__['2']
    if not os.path.isfile(reads2):
        sys.stderr.write('Reads file {} not found.'.format(reads2))
        print
        exit(1)

if not os.path.isfile(reads1):
    sys.stderr.write('Reads file {} not found.'.format(reads1))
    print
    exit(1)

if not os.path.isfile(args.barcodes):
    sys.stderr.write('Barcodes file {} not found.'.format(args.barcodes))
    print
    exit(1)

if args.mapping and not os.path.isfile(args.mapping):
    sys.stderr.write('Barcode mapping file {} not found.'.format(args.mapping))
    print
    exit(1)

pool_maps = {}


def read_line(t, index, code_col):
    try:
        if t[index]:
            pool_maps[t[code_col]] = t[index]
    except:
        pass


def isnum(a):
    try:
        b = int(a)
        return True
    except:
        return False


if args.mapping:
    for l in open(args.mapping, 'rt'):
        l = l.strip()
        if l:
            t = l.split('\t')
            if isnum(t[0]):
                read_line(t, 2, 1)

for bar, map in pool_maps.iteritems():
    print '{}: {}'.format(bar, map)


def read_multiplexed(reads1_file, reads2_file, barcodes_file, pool_maps, progress_start, prefix='temp'):
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    pool_name = reads1_file.split('.')[0]

    def nicename(a):
        return a.replace('#', '').replace('  ', ' ').replace('/', ' ').replace(' ', '_')

    files, f1, f2, fbar = {}, None, None, None
    try:
        barcodes = set(pool_maps.keys())
        print "BARCODES: {}".format(barcodes)

        for barcode in barcodes:
            name = nicename(pool_maps[barcode])
            if reads2_file:
                filename = os.path.join(prefix, '{}_{}_{}_mate1.fq.gz'.format(pool_name, barcode, name))
                files[barcode] = gzip.open(filename, 'wb')

                filename = os.path.join(prefix, '{}_{}_{}_mate2.fq.gz'.format(pool_name, barcode, name))
                files[barcode + '2'] = gzip.open(filename, 'wb')

            else:
                filename = os.path.join(prefix, '{}_{}_{}.fq.gz'.format(pool_name, barcode, name))
                files[barcode] = gzip.open(filename, 'wb')

        if reads2_file:
            files['notmatched'] = gzip.open(
                os.path.join(prefix, '{}_Not_Matched_mate1.fq.gz'.format(pool_name)), 'wb')
            files['badquality'] = gzip.open(
                os.path.join(prefix, '{}_Bad_Quality_mate1.fq.gz'.format(pool_name)), 'wb')
            files['notmatched2'] = gzip.open(
                os.path.join(prefix, '{}_Not_Matched_mate2.fq.gz'.format(pool_name)), 'wb')
            files['badquality2'] = gzip.open(
                os.path.join(prefix, '{}_Bad_Quality_mate2.fq.gz'.format(pool_name)), 'wb')
        else:
            files['notmatched'] = gzip.open(os.path.join(prefix, '{}_Not_Matched.fq.gz'.format(pool_name)), 'wb')
            files['badquality'] = gzip.open(os.path.join(prefix, '{}_Bad_Quality.fq.gz'.format(pool_name)), 'wb')

        filenames = [f.name for f in files.values()]

        p = subprocess.Popen(
            'pigz -dc {} | wc -l'.format(barcodes_file),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        numlines, err = p.communicate()
        if err:
            raise Exception(err)

        numlines = int(numlines)
        id, matched, notmatched, badquality, skipped = 0, 0, 0, 0, 0
        print json.dumps({'proc.progress': progress_start})
        progress = progress_start
        progress_step = (0.9 - progress) / 20.
        progress_span = numlines / 20

        def save_results(matched, notmatched, badquality, skipped, total, progress):
            total = float(total)

            print json.dumps({
                'matched': '{:,} reads ({:.2f} %)'.format(matched, 100 * matched / total),
                'notmatched': '{:,} reads ({:.2f} %)'.format(notmatched, 100 * notmatched / total),
                'badquality': '{:,} reads ({:.2f} %)'.format(badquality, 100 * badquality / total),
                'skipped': '{:,} reads ({:.2f} %)'.format(skipped, 100 * skipped / total),
                'proc.progress': progress
            }, separators=(',', ':'))

        f1 = gzip.GzipFile(reads1_file, 'r')
        r1 = f1.readline()
        fbar = gzip.GzipFile(barcodes_file, 'r')
        rbar = fbar.readline()

        if reads2_file:
            f2 = gzip.GzipFile(reads2_file, 'r')
            r2 = f2.readline()

        while r1:
            id += 1
            r1 = r1.rstrip('\r').rstrip('\n').split('\t')
            s1 = r1[-3].replace('.', 'N')
            rbar = rbar.rstrip('\r').rstrip('\n').split('\t')
            sbar = rbar[-3].replace('.', 'N')[:-1]
            p1, pbar = r1[-1], rbar[-1]

            if reads2_file:
                r2 = r2.rstrip('\r').rstrip('\n').split('\t')
                s2 = r2[-3].replace('.', 'N')
                p2 = r2[-1]
            else:
                r2 = r1
                p2 = p1

            if r1[:7] == r2[:7] == rbar[:7] and p1 == p2 == pbar:
                plusline = '+' + ':'.join(r1[:7]) + ' ' + sbar
                if p1 == '1' and p2 == '1':
                    if sbar in barcodes:
                        files[sbar].write(
                            '@' + str(id) + '\n' + s1 + '\n' + plusline + '\n' + r1[-2] + '\n')
                        if reads2_file:
                            files[sbar + '2'].write(
                                '@' + str(id) + '\n' + s2 + '\n' + plusline + '\n' + r2[-2] + '\n')
                        matched += 1
                    else:
                        files['notmatched'].write(
                            '@' + str(id) + '\n' + s1 + '\n' + plusline + '\n' + r1[-2] + '\n')
                        if reads2_file:
                            files['notmatched2'].write(
                                '@' + str(id) + '\n' + s2 + '\n' + plusline + '\n' + r2[-2] + '\n')
                        notmatched += 1
                else:
                    files['badquality'].write(
                        '@' + str(id) + '\n' + s1 + '\n' + plusline + '\n' + r1[-2] + '\n')
                    if reads2_file:
                        files['badquality2'].write(
                            '@' + str(id) + '\n' + s2 + '\n' + plusline + '\n' + r2[-2] + '\n')
                    badquality += 1
            else:
                print "SKIPPED: {}, p1: {}, p2: {}, pbar: {}".format(id, p1, p2, pbar)
                print "{} ? {} ? {}".format(r1[:7], r2[:7], rbar[:7])
                skipped += 1

            if id % progress_span == 0:
                progress += progress_step
                save_results(matched, notmatched, badquality, skipped, id, progress)

            r1 = f1.readline()
            rbar = fbar.readline()

            if reads2_file:
                r2 = f2.readline()

        save_results(matched, notmatched, badquality, skipped, id, 0.9)

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
        if name.endswith('_mate2'):
            continue

        name2 = name.replace('_mate1', '_mate2')
        d = {
            'status': 'resolving',
            'processor_name': 'import:upload:reads-fastq-paired-end',
            'input': {
                'src1': {
                    'file': os.path.basename(name),
                    'file_temp': os.path.join(os.getcwd(), 'temp', os.path.basename(name))
                },
                'src2': {
                    'file': os.path.basename(name2),
                    'file_temp': os.path.join(os.getcwd(), 'temp', os.path.basename(name2))
                }
            }
        }
    else:
        d = {
            'status': 'resolving',
            'processor_name': 'import:upload:reads-fastq',
            'input': {
                'src': {
                    'file': os.path.basename(name),
                    'file_temp': os.path.join(os.getcwd(), 'temp', os.path.basename(name))
                }
            }
        }

    print 'run {}'.format(json.dumps(d, separators=(',', ':')))

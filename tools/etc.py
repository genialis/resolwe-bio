#!/usr/bin/env python2
import argparse
import gzip
import json
import sys
import numpy as np

from collections import Counter


parser = argparse.ArgumentParser(description='Build expression time course.')
parser.add_argument('expression_files', nargs='+', help='gene expression files')
parser.add_argument('--names', nargs='+', help='expression names')
parser.add_argument('--mean', action='store_true', help='group by time and average')
args = parser.parse_args()

if args.names and len(args.names) != len(args.expression_files):
    raise AttributeError("Names length must match files length.")

files = args.expression_files
progress = 0
series = {'timePoints': [], 'genes': {}}

# Read time points from file names
times = []
if args.names:
    for name in args.names:
        ndx = name.lower().find('hr')
        if ndx > 0:
            try:
                times.append(int(name[ndx + 2:ndx + 4]))
            except:
                try:
                    times.append(int(name[ndx - 2:ndx]))
                except:
                    pass
    if len(times) > 0 and len(times) != len(files):
        raise AttributeError("Found some times but not for all time points.")

if len(times) == 0:
    times = range(len(files))


def is_gzipped(f):
    with open(f, 'rb') as rpkm_file:
        magic = rpkm_file.read(2)

    return magic == '\037\213'


# Go through files
for t, f in sorted(zip(times, files)):
    myopen = gzip.open if is_gzipped(f) else open

    with myopen(f) as rpkm_file:
        rpkm_file.readline()  # skip first line

        # TODO: see above
        series['timePoints'].append(t)

        # Go through lines and parse them
        for line in rpkm_file:
            entry = line.split()
            gene = entry[0]
            value = float(entry[1])

            # Add gene if missing
            if gene not in series['genes']:
                series['genes'][gene] = []

            # Set values from previous experiment to 0 if missing and add value
            series['genes'][gene].extend([0] * (len(series['timePoints']) - 1 - len(series['genes'][gene])))
            series['genes'][gene].append(value)

    # Print progress
    progress += 1.0 / len(files)
    print '{{"proc.progress":{0}}}'.format(progress)

# Add 0s at the end of genes if values are missing
for gene in series['genes']:
    series['genes'][gene].extend([0] * (len(series['timePoints']) - len(series['genes'][gene])))

# Average time points for replicates
if args.mean:
    tp, tp_counts = zip(*sorted(Counter(series['timePoints']).iteritems()))
    for gene in series['genes']:
        s = 0
        tp_averages = []
        for c in tp_counts:
            tp_averages.append(np.mean(series['genes'][gene][s:s + c]))
            s += c

        series['genes'][gene] = tp_averages
    series['timePoints'] = tp

# Print result in json
etcjson = '{"etc":%s}' % json.dumps(series, separators=(',', ':'))
print etcjson
gzip.open('etc.json.gz', 'wb').write(etcjson)

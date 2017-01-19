#!/usr/bin/env python
"""
Insert additional string into filename.

Script solves the problem of inserting additional word into filename between the
actual file name and it's extension. Example use:

    $ extend_name.py -f reads.fasta.gz -e _filtered -v fa fasta
    $ reads_filtered.fasta.gz

"""
from __future__ import absolute_import, division, print_function

import os
import argparse

parser = argparse.ArgumentParser(
    description='Insert additional word between file basename and file extension.')

parser.add_argument('-f', '--filename', help='Original filename (can be a full path).')
parser.add_argument('-a', '--addition', help='Word to be added to the basename.')
parser.add_argument('-e', '--extensions', help='List of possible extensions.', nargs='+')
parser.add_argument('-r', '--remove_extension', help='Remove the extension.', action='store_true')
args = parser.parse_args()

result = []
found = False
for part in os.path.basename(args.filename).split('.')[::-1]:
    if part in args.extensions and not found:
        if not args.remove_extension:
            result = ['.'] + [part] + result
        result = [args.addition] + result
        found = True
    else:
        if not found and args.remove_extension:
            continue
        else:
            result = ['.'] + [part] + result

print(''.join(result[1:]))  # Remove the the first element (['.']) in result

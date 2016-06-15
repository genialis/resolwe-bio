#!/usr/bin/env python2
import argparse
import os
import shutil
import tarfile

parser = argparse.ArgumentParser(description='Decompress and format genome index files.')
parser.add_argument('index', help='Index files (tar.gz bundle)')
parser.add_argument('genome_name', help='Genome base name')
parser.add_argument('output', help='output folder')
parser.add_argument('--bwa', action='store_true', help='Set index file type to BWA')
args = parser.parse_args()

if not tarfile.is_tarfile(args.index):
    raise ValueError('Input file is not a valid .tar.gz file.')

with tarfile.open(args.index) as tar:
    tar.extractall(path='temp_output')

if not os.path.isdir(args.output):
    os.makedirs(args.output)

for root, dirs, files in os.walk('temp_output'):
    for file in files:
        if file.startswith('.'):
            continue
        if args.bwa:
            fname_extention = [file.split('.')[-1]]
        else:
            fname_extention = file.split('.')[1:]
        new_fname = '.'.join([args.genome_name] + fname_extention)
        org_file = os.path.join(root, file)
        new_file = os.path.join(args.output, new_fname)
        shutil.move(org_file, new_file)

shutil.rmtree('temp_output')
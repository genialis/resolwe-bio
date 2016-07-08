#!/usr/bin/env python2
import argparse
import json
import pandas as pd


parser = argparse.ArgumentParser(description="Parse Diff Exp output files")
parser.add_argument('de_file', help="DE analysis output file (.tab).")
args = parser.parse_args()

de_data = pd.read_csv(args.de_file, sep='\t')
de_data.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
columns = {}

for column in de_data:
    columns[column] = list(de_data[column])

data = {'de_data': columns}
print json.dumps(data, separators=(',', ':'))

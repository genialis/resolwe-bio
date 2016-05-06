#!/usr/bin/env python2
import argparse
import json
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Parse Cuffnorm output")
parser.add_argument('input_file', help="genes.fpkm_table file to parse")
parser.add_argument('index', help="column index")
parser.add_argument('-n', action='store_true', help="Normalized input")
args = parser.parse_args()

data = pd.read_csv(args.input_file, sep='\t')
names = [i for i in list(data.columns.values)]
if not args.n:
    expression = pd.concat([data.iloc[:, 0], data.iloc[:, int(args.index) + 1]], axis=1)
    expression.columns = ['Gene', 'Expression']
    expression.to_csv("expression.tab", sep="\t", index=False)
else:
    data['Gene'] = data.index
    expression = pd.concat([data.iloc[:, -1], data.iloc[:, int(args.index)]], axis=1)
    expression.columns = ['Gene', 'Expression']
    expression.to_csv("expression.tab", sep="\t", index=False)

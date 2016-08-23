#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Parse Diff Exp output files."""
from __future__ import absolute_import, division, print_function

import argparse
import json

import numpy as np  # pylint: disable=import-error
import pandas as pd  # pylint: disable=import-error


parser = argparse.ArgumentParser(description="Parse Diff Exp output files")
parser.add_argument('de_file', help="DE analysis output file (.tab).")
args = parser.parse_args()

de_data = pd.read_csv(args.de_file, sep='\t')
de_data.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
de_data.fillna(value=0, inplace=True)
columns = {}

for column in de_data:
    # replace inf values in cuffdiff results file
    if column == 'log2(fold_change)':
        col = np.array(de_data[column])
        col[np.isinf(col)] = 0
        columns[column] = list(col)
    else:
        columns[column] = list(de_data[column])

data = {'de_data': columns}
print(json.dumps(data, separators=(',', ':'), allow_nan=False))

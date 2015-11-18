#!/usr/bin/env python2
import argparse
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

parser = argparse.ArgumentParser(description="Parse Kallisto output")
parser.add_argument('input_file', help="Kallisto output (.tsv) file.")
args = parser.parse_args()

data = pd.read_csv(args.input_file, sep='\t')

tpm = data[['target_id', 'tpm']]
tpm.columns = ['Gene', 'Expression']

rc = data[['target_id', 'est_counts']]
rc['est_counts'] = rc['est_counts'].round()
rc.columns = ['Gene', 'Expression']

tpm.to_csv("abundance_tpm.tab", sep="\t", index=False)
rc.to_csv("abundance_rc.tab", sep="\t", float_format='%.f', index=False)

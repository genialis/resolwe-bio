#!/usr/bin/env python2
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Parse Cuffnorm output")
parser.add_argument('input_file', help="isoforms.fpkm_tracking file to parse")
args = parser.parse_args()

data = pd.read_csv(args.input_file, sep='\t')
names = ["class_code", "nearest_ref_id"] + [i for i in list(data.columns.values) if '_FPKM' in i]
expression = pd.concat([data[i] for i in names], axis=1)
out_data = expression[expression.class_code == "="]
out_data = out_data.drop("class_code", 1)

out_data.to_csv("expression_set.tsv", sep="\t", index=False)

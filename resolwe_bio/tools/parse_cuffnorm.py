#!/usr/bin/env python3
# XXX: Refactor to a comand line tool and remove pylint disable
"""Parse Cuffnorm output."""
import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Parse Cuffnorm output")
parser.add_argument("input_file", help="genes.fpkm_table file to parse")
parser.add_argument("sample_name", help="Sample name")
parser.add_argument("-n", action="store_true", help="Normalized input")
args = parser.parse_args()

data = pd.read_csv(args.input_file, sep="\t")
names = [value.rsplit("_Group", 1)[0] for value in list(data.columns.values)]
index = [i for i, name in enumerate(names) if args.sample_name == name][0]
if not args.n:
    expression = pd.concat([data.iloc[:, 0], data.iloc[:, int(index) + 1]], axis=1)
    expression.columns = ["Gene", "Expression"]
    expression.to_csv("expression.tab", sep="\t", index=False)
else:
    data["Gene"] = data.index
    expression = pd.concat([data.iloc[:, -1], data.iloc[:, int(index)]], axis=1)
    expression.columns = ["Gene", "Expression"]
    expression.to_csv("expression.tab", sep="\t", index=False)

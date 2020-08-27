#!/usr/bin/env python3
"""Import expression time course."""

import csv
import gzip
import json
import sys

import xlrd

file_name = sys.argv[1]


def import_excel(file_name):
    """Import expression time course from Excel."""
    genes = {}
    workbook = xlrd.open_workbook(file_name)
    worksheet = workbook.sheets()[0]
    num_rows = worksheet.nrows
    num_cells = worksheet.ncols
    times = [worksheet.cell_value(0, x) for x in range(1, num_cells)]
    for x in range(1, num_rows):
        genes[str(worksheet.cell_value(x, 0))] = [
            worksheet.cell_value(x, y) for y in range(1, num_cells)
        ]
    return times, genes


def import_table(file_name):
    """Import expression time course from tab separated file."""
    genes = {}
    with open(file_name, "rb") as csvfile:
        my_reader = csv.reader(csvfile, delimiter="\t")
        times = [int(x) for x in my_reader.next()[1:]]
        for x in my_reader:
            genes[str(x[0])] = x[1:]
    return times, genes


if file_name[-4:] == ".xls" or file_name[-5:] == ".xlsx":
    times, genes = import_excel(file_name)

else:
    times, genes = import_table(file_name)

etcjson = '{"etc":%s}' % json.dumps(
    {"genes": genes, "timePoints": times}, separators=(",", ":")
)
print(etcjson)
zipfile = gzip.GzipFile(
    filename="",
    mode="wb",
    fileobj=open("etc.json.gz", "wb"),
    mtime=0,
)
zipfile.write(etcjson.encode("utf-8"))

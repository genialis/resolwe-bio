#!/usr/bin/env python2
import sys
import xlrd
import csv
import gzip
import json

file_name = sys.argv[1]


def import_excel(file_name):
    genes = {}
    workbook = xlrd.open_workbook(file_name)
    worksheet = workbook.sheets()[0]
    num_rows = worksheet.nrows
    num_cells = worksheet.ncols
    times = [worksheet.cell_value(0, x) for x in range(1, num_cells)]
    for x in range(1, num_rows):
        genes[str(worksheet.cell_value(x, 0))] = [worksheet.cell_value(x, y) for y in range(1, num_cells)]
    return times, genes


def import_table(file_name):
    genes = {}
    with open(file_name, 'rb') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        times = [int(x) for x in my_reader.next()[1:]]
        for x in my_reader:
            genes[str(x[0])] = x[1:]
    return times, genes
if file_name[-4:] == ".xls" or file_name[-5:] == ".xlsx":
    times, genes = import_excel(file_name)

else:
    times, genes = import_table(file_name)

etcjson = '{"etc":%s}' % json.dumps({'genes': genes, 'timePoints': times}, separators=(',', ':'))
print etcjson
gzip.open('etc.json.gz', 'wb').write(etcjson)

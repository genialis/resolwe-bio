#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,redefined-outer-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Import gene information."""
from __future__ import absolute_import, division, print_function

import csv
import json
import os
import sys

import xlrd  # pylint: disable=import-error


file_name = sys.argv[1]
if not os.path.isfile(file_name):
    raise ValueError("File {} does not exist".format(file_name))


def import_excel(file_name):
    """Import gene information from Excel."""
    meta = {'data': [], "header": []}
    workbook = xlrd.open_workbook(file_name)
    worksheet = workbook.sheets()[0]
    num_rows = worksheet.nrows
    num_cols = worksheet.ncols

    for x in range(0, num_cols):
        meta["header"].append(str(worksheet.cell_value(0, x)))

    data = {}

    for i, l in enumerate(meta["header"]):
        data[l] = i

    for x in range(1, num_rows):
        out = {}

        for label in data:
            out[label] = str(worksheet.cell_value(x, data[label]))

        meta['data'].append(out)
    return meta


def import_table(file_name):
    """Import gene information from tab separated file."""
    meta = {'data': [], 'header': []}
    with open(file_name, 'rU') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        header = next(my_reader, None)
        meta["header"] = header
        data = {}

        for i, l in enumerate(header):
            data[l] = i

        for vals in my_reader:
            out = {}

            for label in data:
                out[label] = vals[data[label]]

            meta['data'].append(out)
    return meta


if file_name[-4:] == ".xls" or file_name[-5:] == ".xlsx":
    meta = import_excel(file_name)
else:
    meta = import_table(file_name)


print(json.dumps({'info': meta}, separators=(',', ':')))

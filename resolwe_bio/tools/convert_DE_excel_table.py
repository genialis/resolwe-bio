#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Conver differential expression to Excel table."""
from __future__ import absolute_import, division, print_function

import sys
import xlrd  # pylint: disable=import-error
file_name = sys.argv[1]


workbook = xlrd.open_workbook(file_name)
worksheet = workbook.sheets()[0]

num_rows = worksheet.nrows
curr_row = 0

print("\t".join([i.value for i in worksheet.row(0)]))

for curr_row in range(1, worksheet.nrows):
    print("\t".join([str(i.value) for i in worksheet.row(curr_row)]))

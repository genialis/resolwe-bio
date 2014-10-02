import sys
import xlrd
import csv
import json
import os

file_name = sys.argv[1]
if not os.path.isfile(file_name):
    raise ValueError("File {} does not exist".format(file_name.gene_information_file))


def import_excel(file_name):
    meta = {'geneinfo': []}
    workbook = xlrd.open_workbook(file_name)
    worksheet = workbook.sheets()[0]
    num_rows = worksheet.nrows
    for x in range(1, num_rows):
        meta['geneinfo'].append({
            "id": str(worksheet.cell_value(x, 0)),
            "name": str(worksheet.cell_value(x, 1)),
            "synonyms": str(worksheet.cell_value(x, 2)),
            "description": str(worksheet.cell_value(x, 3))})
    return meta


def import_table(file_name):
    meta = {'geneinfo': []}
    with open(file_name, 'rU') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        next(my_reader, None)
        for vals in my_reader:
            meta['geneinfo'].append({"id": vals[0], "name": vals[1],
                                     "synonyms": vals[2],
                                     "description": vals[3]})
    return meta


if file_name[-4:] == ".xls" or file_name[-5:] == ".xlsx":
    meta = import_excel(file_name)
else:
    meta = import_table(file_name)


print json.dumps({'info': meta}, separators=(',', ':'))

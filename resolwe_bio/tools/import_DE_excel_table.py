import sys
import xlrd
file_name = sys.argv[1]


def check_headers(header):
    """
    Checks header if arguments are valid: full string match or start string match
    :param header:
    :return: True if all test pass
    """
    valid_titles = ['Gene', 'Case_Counts_<identifier>', 'Control_Counts_<identifier>',
                    'Lik.DE', 'FDR.DE', 'Lik.NDE', 'FDR.NDE', 'Case_RPKUM_Median', 'Control_RPKUM_Median']
    tests = [
        header[0].value == 'Gene',
        header[1].value.startswith('Case_Counts_'),
        header[2].value.startswith('Control_Counts_'),
        header[3].value == 'Lik.DE',
        header[4].value == 'FDR.DE',
        header[5].value == 'Lik.NDE',
        header[6].value == 'FDR.NDE',
        header[7].value == 'Case_RPKUM_Median',
        header[8].value == 'Control_RPKUM_Median']

    if all(tests):
        return True
    invalid = tests.index(False)
    sys.stderr.write('##########\nERROR: invalid header at column number {}\nERROR: Index should be'
                     ' {}\n##########\n'.format(invalid + 1, valid_titles[invalid]))
    sys.exit('exiting with error')


def check_values(row, nu):
    """
    Checks if all rows are valid (column 1: not empty, column 2-9: float )
    :param row:
    :return:
    """
    if all([x.ctype == 2 for x in row[1:10]]) and row[0].ctype != 0:
        return True
    sys.stderr.write('##########\nERROR: invalid value at row number {}\nERROR: First field should not'
                     'be empty and other fileds should contain numbers\n##########\n'.format(nu + 1))
    sys.exit("exiting with error")

workbook = xlrd.open_workbook(file_name)
worksheet = workbook.sheets()[0]

num_rows = worksheet.nrows
curr_row = 0
check_headers(worksheet.row(0))

print "\t".join([i.value for i in worksheet.row(0)])


for curr_row in range(1, worksheet.nrows):
    check_values(worksheet.row(curr_row), curr_row)
    print "\t".join([str(i.value) for i in worksheet.row(curr_row)])

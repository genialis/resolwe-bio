import sys
import csv
file_name = sys.argv[1]


def check_headers(header):
    """
    Checks header if arguments are valid: full string match or start string match
    :param header:
    :return: True if all test pass
    """
    valid_titles = ['Gene', 'Case_Counts_<identifier>', 'Control_Counts_<identifier>', 'Lik.DE',
                    'FDR.DE', 'Lik.NDE', 'FDR.NDE', 'Case_RPKUM_Median', 'Control_RPKUM_Median']
    tests = [
        header[0] == 'Gene',
        header[1].startswith('Case_Counts_'),
        header[2].startswith('Control_Counts_'),
        header[3] == 'Lik.DE',
        header[4] == 'FDR.DE',
        header[5] == 'Lik.NDE',
        header[6] == 'FDR.NDE',
        header[7] == 'Case_RPKUM_Median',
        header[8] == 'Control_RPKUM_Median']

    if all(tests):
        return True
    invalid = tests.index(False)
    sys.stderr.write('##########\nERROR: invalid header at column number {}\nERROR: '
                     'Index should be "{}"\n##########\n'.format(invalid + 1, valid_titles[invalid]))
    sys.exit('exiting with error')


def check_values(row, nu):
    """
    Checks if all rows are valid (column 1: not empty, column 2-9: float )
    :param row:
    :return: TRUE if all tests are passed, else ERROR
    """
    try:
        assert len(row) >= 9
        assert row[0]
        [float(x) for x in row[1:9]]
    except:
        sys.stderr.write('##########\nERROR: invalid value at row number {}\nERROR: First field'
                         'should not be empty and other fileds should contain numbers\n##########\n'.format(nu + 1))
        sys.exit('exiting with error')


curr_row = 0
with open(file_name) as tsv:
    r_hander = csv.reader(tsv, delimiter='\t')
    check_headers(r_hander.next())
    for line in r_hander:
        curr_row += 1
        check_values(line, curr_row)

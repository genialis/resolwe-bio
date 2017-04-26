#!/usr/bin/env python3
"""Report iClip sample QC status."""
import argparse
import json


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Report iClip sample QC status.")
    parser.add_argument('reads', help="Number of reads in the alignment (BAM) file.")
    parser.add_argument('cdna', help="Number of cDNAs.")
    return parser.parse_args()


def report_qc(cdna, reads):
    """Report sample QC."""
    ratio = reads / cdna
    status, message = 'FAIL', ''

    if cdna > 500000:
        status = 'PASS'
        message = 'Sample contains more then 500,000 aligned reads.'
    elif cdna > 200000 and cdna <= 500000 and ratio < 3:
        status = 'WARNING'
        message = 'Sample contains less then 500,000 aligned reads, but the reads/cdna ratio is below 3.'
    elif cdna <= 200000:
        message = 'Low number of aligned reads, less then 200,000.'
    else:
        message = 'The number of reads is between 200,000 and 500,000, but the reads/cdna ratio is above 3.'

    return status, message


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    status, message = report_qc(int(args.cdna), int(args.reads))
    print(json.dumps({'qc': {'status': status, 'message': message}}, separators=(',', ':')))


if __name__ == "__main__":
    main()

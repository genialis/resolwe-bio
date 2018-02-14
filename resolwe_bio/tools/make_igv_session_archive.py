#!/usr/bin/env python3
"""Create IGV session for archive."""

import argparse
import os

from lxml import etree

parser = argparse.ArgumentParser(description="Create igv session.")

parser.add_argument('-f', '--input_file', required=True, help="File with paths to files for IGV")

args = parser.parse_args()


def make_xml_tree(input_file):
    """Make xml tree for IGV session."""
    global_ = etree.Element(  # pylint: disable=no-member
        'Global',
        version='3',
    )

    resources = etree.SubElement(global_, 'Resources')  # pylint: disable=no-member

    with open(input_file, 'r') as tfile:
        for line in tfile:
            line = line.rstrip()
            etree.SubElement(  # pylint: disable=no-member
                resources,
                'Resource',
                name=os.path.basename(line),
                path=os.path.join('..', line),
            )

    doc = etree.tostring(  # pylint: disable=no-member
        global_,
        pretty_print=True,
        xml_declaration=True,
        encoding='UTF-8',
    )
    return doc


def write_xml_file(input_file, doc):
    """Compose a name and write output file."""
    with open(input_file.replace('temp_igv.txt', 'igv.xml'), 'wb') as f:
        f.write(doc)


def main():
    """Invoke when run directly as a program."""
    doc = make_xml_tree(args.input_file)
    write_xml_file(args.input_file, doc)


if __name__ == "__main__":
    main()

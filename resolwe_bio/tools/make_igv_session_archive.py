#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name,import-error
"""Create IGV session for archive."""

import argparse
from lxml import etree
import os

parser = argparse.ArgumentParser(description='Create igv session')

parser.add_argument("-f", "--files", nargs='+', required=True, help="List of files")
parser.add_argument("-g", "--genome", type=str, required=True, help="Genome")

args = parser.parse_args()

Global = etree.Element(
    'Global',
    genome=args.genome,
    version='8'
)

Resources = etree.SubElement(Global, 'Resources')

for element in args.files:
    element_name = os.path.basename(element)
    element_dir = os.path.dirname(os.path.dirname(element))
    element_path = os.path.join('..', element)

    if element_dir == 'None':
        element_subdir = os.path.basename(os.path.dirname(element))
        element_path = os.path.join('..', 'other_data', element_subdir, element_name)

    Resource = etree.SubElement(Resources, 'Resource', name=element_name, path=element_path)

doc = etree.tostring(Global, pretty_print=True, xml_declaration=True, encoding='UTF-8')

with open('igv_session.xml', 'wb') as f:
    f.write(doc)

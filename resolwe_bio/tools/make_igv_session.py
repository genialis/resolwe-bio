#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name,import-error
"""Create IGV session."""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
from lxml import etree

parser = argparse.ArgumentParser(description='Create igv session')

parser.add_argument("-bam", "--bam_files", nargs='+', required=True, help="List of paths to bam files")
parser.add_argument("-g", "--genomeid", required=True, help="Genome id")
parser.add_argument("-l", "--bam_locus", required=True, help="Locus where to open IGV")
args = parser.parse_args()

Global = etree.Element(
    'Global',
    genome=args.genomeid,
    hasGeneTrack='true',
    hasSequenceTrack='true',
    locus=args.bam_locus,
    nextAutoscaleGroup='3',
    version='8'
)

Resources = etree.SubElement(Global, 'Resources')

for bam_path in args.bam_files:
    Resources.append(etree.Element('Resource', path=bam_path))

doc = etree.tostring(Global, pretty_print=True, xml_declaration=True, encoding='UTF-8')

with open('xmltree.xml', 'w') as f:
    f.write(doc)

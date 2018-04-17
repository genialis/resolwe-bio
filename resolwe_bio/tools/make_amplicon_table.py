#!/usr/bin/env python2
"""Create amplicon variant table."""
from __future__ import absolute_import, division, print_function

from collections import defaultdict

import argparse
import json

import HTSeq


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Create amplicon variant table.")
    parser.add_argument('master_file', help="Amplicon master file.")
    parser.add_argument('amplicon_coverage', help='Per-amplicon mean coverage.')
    parser.add_argument('-v', '--variants', nargs='+', required=True, help='snpEff output file')
    parser.add_argument('-a', '--list_all', action='store_true', help='Output all amplicons.')
    return parser.parse_args()


def _tsv_to_list(table_file, has_header=False, delimiter='\t', pick_columns=None):
    """Parse *.tsv file into list."""
    table_data = []
    header = None
    with open(table_file, 'r') as tfile:
        if has_header:
            header = next(tfile).strip().split(delimiter)
            common_columns = [x for x in pick_columns if x in header]
            if pick_columns:
                # Find indexes of selected columns
                temp_header = {col: i for i, col in enumerate(header)}
                pick_indexes = [temp_header[col] for col in common_columns]
                header = common_columns
        for line in tfile:
            line_content = line.strip().split(delimiter)
            if pick_columns:
                # Select only specific columns and order them as in pick_columns:
                temp_line_content = {i: col for i, col in enumerate(line_content)}
                line_content = [temp_line_content[i] for i in pick_indexes]
            table_data.append(line_content)
    return table_data, header


def snp_url(snpid):
    """Create url SNP ID."""
    if snpid.startswith('rs'):
        url = 'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
        pass
    elif snpid.startswith('COSM'):
        url = 'http://cancer.sanger.ac.uk/cosmic/mutation/overview?genome=37\&id={}'.format(snpid.lstrip('COSM'))
    elif snpid.startswith('COSN'):
        url = 'http://cancer.sanger.ac.uk/cosmic/ncv/overview?genome=37\&id={}'.format(snpid.lstrip('COSN'))
    else:
        return snpid
    return url


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    amplicons = {}
    amp_cov = {}
    variants_temp = defaultdict(list)
    variants = []

    with open(args.master_file) as master_file:
        for line in master_file:
            amplicon = line.strip().split('\t')

            amplicons[amplicon[3]] = HTSeq.GenomicInterval(
                amplicon[0],
                int(amplicon[1]),
                int(amplicon[2])
            )

    with open(args.amplicon_coverage) as ac:
        for line in ac:
            cov = line.strip().split('\t')
            amp_cov[cov[0]] = float(cov[1])

    for annot_var in args.variants:
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP', 'DP4', 'GEN[0].AD', 'SB', 'FS', 'EFF[*].GENE', 'ID']
        data, header = _tsv_to_list(annot_var, has_header=True, pick_columns=columns)

        for var in data:
            variant = dict(zip(header, var))
            gb_pos = '{}:{}'.format(variant['CHROM'], int(variant['POS']))
            var_or_indel = '{}>{}'.format(variant['REF'], variant['ALT'])

            variant_pos = HTSeq.GenomicPosition(variant['CHROM'], int(variant['POS']))
            for amp in amplicons:
                if variant_pos.is_contained_in(amplicons[amp]):
                    cov = amp_cov[amp]
                    break

            if 'DP4' in variant:
                process = 'LoFreq'
                other_fields = 'DP4={};SB={}'.format(variant['DP4'], variant['SB'])
            else:
                process = 'GATK'
                other_fields = 'AD:{};FS:{}'.format(variant['GEN[0].AD'], variant['FS'])

            if variant['ID'].strip(';') == '':
                var_links = []
            else:
                var_links = [(var_id, snp_url(var_id)) for var_id in variant['ID'].split(';')]

            url = 'http://www.ncbi.nlm.nih.gov/gene/?term='
            feature_links = [(variant['EFF[*].GENE'], '{}{}'.format(url, variant['EFF[*].GENE']))]

            variants_temp[amp].append({
                'pos': gb_pos,
                'columns': [
                    amp, str(cov), feature_links, variant['CHROM'], variant['POS'], var_or_indel,
                    variant['AF'], variant['DP'], other_fields, process, var_links
                ]
            })

    for amp in amplicons:
        if amp in variants_temp:
            for variant in variants_temp[amp]:
                variants.append({
                    'pos': variant['pos'],
                    'columns': variant['columns']
                })
        elif amp not in variants_temp and args.list_all:
            variants.append({
                'pos': '{}:{}..{}'.format(amplicons[amp].chrom, amplicons[amp].start, amplicons[amp].end),
                'columns': [amp, str(amp_cov[amp]), '', '', '', '', '', '', '', '', '']
            })

    column_types = [
        'value',
        'value',
        'urls',
        'value',
        'value',
        'value',
        'value',
        'value',
        'delimited',
        'value',
        'urls'
    ]

    nice_headers = [
        'A. ID',
        'A. Cov.',
        'Gene',
        'CHR',
        'V. POS',
        'REF>ALT',
        'AF',
        'DoC',
        'Other',
        'Tool',
        'Links'
    ]

    header_labels = [
        'Amplicon ID',
        'Amplicon Coverage',
        'Gene symbol',
        'Chromosome',
        'Variant position',
        'Identified variant/indel',
        'Allele frequency',
        'Depth of coverage (filtered depth)',
        'Other variant stats',
        'Tool',
        'Links'
    ]

    amp_table = {
        'headers': nice_headers,
        'labels': header_labels,
        'data': variants,
        'column_types': column_types
    }

    with open('amplicon_table.json', 'w') as f:
        json.dump(amp_table, f)


if __name__ == "__main__":
    main()

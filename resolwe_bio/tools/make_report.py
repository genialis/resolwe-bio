#!/usr/bin/env python2
"""Generate amplicon report."""
from __future__ import division

import os
import argparse
import subprocess


DECIMALS = 2  # Decimal precision when presenting results:

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.")
parser.add_argument('--panel', help="Panel name")
parser.add_argument('--covmetrics', help="Coverge metrics")
parser.add_argument('--cov', help="Amplicon coverage")
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.")
parser.add_argument('--vcf', help="File with VCF data.", nargs='+')
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")


def _escape_latex(string):
    """Format normal string to be LaTeX compliant."""
    string = string.replace('\\', '\\\\')
    string = string.replace('&', '\&')
    string = string.replace('%', '\%')
    string = string.replace('$', '\$')
    string = string.replace('#', '\#')
    string = string.replace('_', '\_')
    string = string.replace('{', '\{')
    string = string.replace('}', '\}')
    string = string.replace('~', '\textasciitilde ')
    string = string.replace('^', '\textasciicircum ')
    return string


def _remove_underscore(string):
    """Remove underscore from the LaTeX compliant string and replaces it with white space."""
    return string.replace('\_', ' ')


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


def cut_to_pieces(string, piece_size):
    """Cut long string into smaller pieces."""
    pieces = []
    while len(string) > piece_size:
        pieces.append(string[:piece_size])
        string = string[piece_size:]
    pieces.append(string)
    return ' '.join(pieces)


def list_to_tex_table(data, header=None, caption=None, long_columns=False):
    """Make a TeX table from python list."""
    lines = []
    column_alingnment = ['l' for h in header]
    if long_columns:
        for col_index in long_columns:
            column_alingnment[col_index] = 'L'

    if caption:
        lines.append('\\captionof{{table}}{{{}}}'.format(caption))

    lines.append('\\noindent')
    lines.append('\\begin{{longtable}}[l] {{ {} }}'.format('  '.join(column_alingnment)))

    if header:
        lines.append('\\rowcolor{darkblue1}')
        lines.append('\\leavevmode\\color{white}\\textbf{' +
                     '}& \\leavevmode\\color{white}\\textbf{'.join(header) + '} \\\\')

    for line in data:
        if long_columns:
            for col_index in long_columns:
                # If hyperlink line, don't do a thing. Otherwise, insert spaces, so that wrapping can happen:
                new_val = line[col_index] if '\href' in line[col_index] else cut_to_pieces(line[col_index], 8)
                line[col_index] = '\\multicolumn{{1}}{{m{{2.3cm}}}}{{{}}}'.format(new_val)

        lines.append(' & '.join(line) + ' \\\\')

    lines.append('\\end{longtable}')
    lines.append('{\n\\addtocounter{table}{-1}}')
    return '\n'.join(lines)


def snp_href(snpid):
    """Create LaTeX hyperlink for given SNP ID."""
    if snpid.startswith('rs'):
        url = 'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
        pass
    elif snpid.startswith('COSM'):
        url = 'http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id={}'.format(snpid.lstrip('COSM'))
    elif snpid.startswith('COSN'):
        url = 'http://cancer.sanger.ac.uk/cosmic/ncv/overview?id={}'.format(snpid.lstrip('COSN'))
    else:
        return snpid
    return '\\href{{{}}}{{{}}}'.format(url, snpid)


def gene_href(gene_name):
    """Create LaTeX hyperlink for given GENE ID."""
    url = 'http://www.ncbi.nlm.nih.gov/gene/?term={}'.format(gene_name)
    return '\\href{{{}}}{{{}}}'.format(url, gene_name)


def parse_target_pcr_metrics(metrics_report):
    """Parse CollectTargetedPcrMetrics report file."""
    with open(metrics_report) as pcr_metrics:
        labels = []

        for line in pcr_metrics:
            if line.startswith('#'):
                continue
            if len(labels) == 0:
                if line.startswith('CUSTOM_AMPLICON_SET'):
                    labels = line.strip().split('\t')
            else:
                values = line.strip().split('\t')
                break

        return dict(zip(labels, values))


def vcf_table_name(vcf_file_name):
    """Format VCF table caption."""
    if 'gatkhc.finalvars' in vcf_file_name.lower():
        return 'GATK HaplotypeCaller variant calls'
    elif 'lf.finalvars' in vcf_file_name.lower():
        return 'Lofreq variant calls'
    else:
        return os.path.basename(vcf_file_name)


if __name__ == '__main__':
    args = parser.parse_args()

    # Open template ind fill it with data:
    with open(args.template, 'r') as template_in, open('report.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)
        template = template.replace('{#SAMPLE_NAME#}', _escape_latex(args.sample))
        template = template.replace('{#PANEL#}', _remove_underscore(_escape_latex(args.panel)))

        # get coverage stats
        with open(args.covmetrics) as covmetrics:
            cov_data = covmetrics.readline().strip().split('\t')
            mean_coverage = float(cov_data[0])
            mean20 = float(cov_data[1])
            cov_unif = float(cov_data[2])

        # Produce results, metrics:
        metrics = parse_target_pcr_metrics(args.metrics)

        pct_aligned_reads = str('{0:g}'.format(round(float(metrics['PCT_PF_UQ_READS_ALIGNED']) * 100, DECIMALS)))
        pct_amplified_bases = str('{0:g}'.format(round(float(metrics['PCT_AMPLIFIED_BASES']) * 100, DECIMALS)))

        template = template.replace('{#TOTAL_READS#}', metrics['TOTAL_READS'])
        template = template.replace('{#ALIGNED_READS#}', metrics['PF_UQ_READS_ALIGNED'])
        template = template.replace('{#PCT_ALIGNED_READS#}', pct_aligned_reads)
        template = template.replace('{#ALIGN_BASES_ONTARGET#}', pct_amplified_bases)
        template = template.replace('{#COV_MEAN#}', str(round(mean_coverage, DECIMALS)))
        template = template.replace('{#COV_20#}', str(round(mean20, DECIMALS)))
        template = template.replace('{#COV_UNI#}', str(round(cov_unif, DECIMALS)))

        # Parse *.cov file:
        cov_list, _ = _tsv_to_list(args.cov)
        amp_numb = list(range(len(cov_list)))
        covered_amplicons = len([1 for line in cov_list if float(line[8]) >= 1])

        # Produce LaTeX table for amplicons with insufficient coverage:
        not_covered_amplicons = [
            (_escape_latex(line[4]), '{:.1f}'.format(float(line[8]) * 100)) for line in cov_list if float(line[8]) < 1]
        switch = 0  # switch used for adding the sentence about coverage

        if not not_covered_amplicons:
            switch = 1
            not_covered_amplicons = [['/', '/']]
        table_text = list_to_tex_table(not_covered_amplicons, header=['Amplicon name', '\% covered'],
                                       caption='Amplicons with coverage $<$ 100\%')

        template = template.replace('{#ALL_AMPLICONS#}', str(len(cov_list)))
        template = template.replace('{#COVERED_AMPLICONS#}', str(covered_amplicons))
        template = template.replace(
            '{#PCT_COV#}', str('{0:g}'.format(round(float(covered_amplicons) / len(cov_list) * 100, DECIMALS))))
        if switch:
            message = 'All amplicons are completely covered with short reads.'
            sentence_text = '\\medskip \\noindent \n {{\\boldfont ' + message + '}} \n\n \\lightfont \n'
            template = template.replace('{#BAD_AMPLICON_TABLE#}', sentence_text + table_text)
        else:
            template = template.replace('{#BAD_AMPLICON_TABLE#}', table_text)

        # Make VCF tables:
        table_text = ''
        for vcf_file in args.vcf:
            header = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP', 'DP4', 'GEN[0].AD', 'SB', 'FS', 'EFF[*].GENE', 'ID']
            vcf_table, common_columns = _tsv_to_list(vcf_file, has_header=True, pick_columns=header)

            # Escape user inputs:
            common_columns = [_escape_latex(name) for name in common_columns]
            caption = _escape_latex(vcf_table_name(vcf_file))
            vcf_table = [[_escape_latex(value) for value in line] for line in vcf_table]

            # Insert space between SNP ID's and create hypelinks:
            vcf_table = [line[:-1] + [' '.join(map(snp_href, line[-1].split(';')))] for line in vcf_table]

            # Create gene hypelinks:
            vcf_table = [line[:-2] + [gene_href(line[-2])] + [line[-1]] for line in vcf_table]

            table_text += list_to_tex_table(vcf_table, header=common_columns, caption=caption, long_columns=[2, 3, -1])
            table_text += '\n\\newpage\n'
        template = template.replace('{#VCF_TABLES#}', table_text)

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'report.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    print("Done.")

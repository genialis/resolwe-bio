#!/usr/bin/env python3
"""Generate amplicon report."""
import argparse
import csv
import os
import re
import subprocess

DECIMALS = 2  # Decimal precision when presenting results:

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.")
parser.add_argument('--panel', help="Panel name")
parser.add_argument('--covmetrics', help="Coverge metrics")
parser.add_argument('--cov', help="Amplicon coverage")
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.")
parser.add_argument('--annot_vars', help="Annotated variant files.", nargs='+')
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")
parser.add_argument('--afthreshold', help="Allele Frequency lower threshold.")


def escape_latex(string):
    """
    Get string, where reserved LaTeX charachters are escaped.

    For more info regarding LaTeX reserved charachters read this:
    https://en.wikibooks.org/wiki/LaTeX/Basics#Reserved_Characters
    """
    string = string.replace('\\', '\\\\')
    string = string.replace('&', r'\&')
    string = string.replace('%', r'\%')
    string = string.replace('$', r'\$')
    string = string.replace('#', r'\#')
    string = string.replace('_', r'\_')
    string = string.replace('{', r'\{')
    string = string.replace('}', r'\}')
    string = string.replace('~', r'\textasciitilde ')
    string = string.replace('^', r'\textasciicircum ')
    return string


def format_float(value, decimals=DECIMALS, to_percent=False):
    """Round string representation of float to ``decimals`` decimal places.

    If ``to_percent``==True, also multiply by 100 and add latex percent sign.
    """
    value = float(value) * 100 if to_percent else float(value)
    out = '{0:g}'.format(round(value, decimals))
    if to_percent:
        out += '\\%'
    return out


def parse_covmetrics(covmetrics_file, stats):
    """
    Parse covmetrics report and write content to ``stats``.

    Covmetrics report is a one-line file with three colunmns::

        mean_coverage   mean_coverage*0.2   coverage_uniformity

    """
    with open(covmetrics_file) as handle:
        cov_data = handle.readline().strip().split('\t')
        stats['mean_coverage'] = float(cov_data[0])
        stats['mean_coverage_20'] = float(cov_data[1])
        stats['coverage_uniformity'] = float(cov_data[2]) / 100  # Convert from percent to ratio


def parse_target_pcr_metrics(metrics_report, stats):
    """
    Parse Picard PcrMetrics report file and save relevant data to ``stats``.

    The are only two relevant lines: the one starting with CUSTOM_AMPLICON_SET
    (containing labels) and the next one (containing coresponding values).
    """
    with open(metrics_report) as pcr_metrics:
        for line in pcr_metrics:
            if line.startswith('#'):
                continue
            if line.startswith('CUSTOM_AMPLICON_SET'):
                labels = line.strip().split('\t')
                values = next(pcr_metrics).strip().split('\t')

        for lab, val in zip(labels, values):
            stats[lab] = val


def cut_to_pieces(string, piece_size):
    """Cut long string into smaller pieces."""
    assert isinstance(string, str)
    string_size = len(string)
    pieces = [string[i: i + piece_size] for i in range(0, string_size, piece_size)]
    return ' '.join(pieces)


def list_to_tex_table(data, header=None, caption=None, long_columns=False, uncut_columns=[]):
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
        header_formatted = ['\\leavevmode\\color{{white}}\\textbf{{{}}}'.format(hdr) for hdr in header]
        lines.append('\\rowcolor{darkblue1}')
        lines.append(' & '.join(header_formatted) + '\\\\')
        lines.append('\\endhead')

    for line in data:
        if long_columns:
            for col_index in long_columns:
                if (r'\href' in line[col_index] or col_index in uncut_columns):
                    # If hyperlink line or `uncut_columns`, don't do a thing.
                    new_val = line[col_index]
                else:
                    # Otherwise, insert spaces, so that wrapping can happen:
                    new_val = cut_to_pieces(line[col_index], 8)
                line[col_index] = '\\multicolumn{{1}}{{m{{2.3cm}}}}{{{}}}'.format(new_val)

        lines.append(' & '.join(line) + ' \\\\')

    lines.append('\\end{longtable}')
    lines.append('{\n\\addtocounter{table}{-1}}')
    return '\n'.join(lines)


def snp_href(snpid):
    """Create LaTeX hyperlink for given SNP ID."""
    if snpid.startswith('rs'):
        url = r'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
    elif snpid.startswith('COSM'):
        url = r'http://cancer.sanger.ac.uk/cosmic/mutation/overview?genome=37\&id={}'.format(snpid.lstrip('COSM'))
    elif snpid.startswith('COSN'):
        url = r'http://cancer.sanger.ac.uk/cosmic/ncv/overview?genome=37\&id={}'.format(snpid.lstrip('COSN'))
    else:
        return snpid
    return '\\href{{{}}}{{{}}}'.format(url, escape_latex(snpid))


def gene_href(gene_name):
    """Create LaTeX hyperlink for given GENE ID."""
    url = 'http://www.ncbi.nlm.nih.gov/gene/?term={}'.format(gene_name)
    return '\\href{{{}}}{{{}}}'.format(url, escape_latex(gene_name))


def format_aa_change(aa_list):
    """Create Amino Acid Change information."""
    output = set()
    if aa_list:
        for aa in aa_list:
            match_obj = re.match(r'p\.([A-Za-z]*)[0-9]*([A-Za-z]*)', aa)
            if match_obj and match_obj.group(1) == match_obj.group(2):
                output.add('Synon')
            else:
                output.add(escape_latex(aa))
    return output


def make_variant_table(variant_file, header, af_threshold):
    """Make variant table."""
    var_table = []
    with open(variant_file, 'r') as vfile:
        for row_raw in csv.DictReader(vfile, delimiter='\t'):
            # Reorder columns according to header
            row = [row_raw[column_name] for column_name in header]

            row[-3] = ' '.join(map(gene_href, row[-3].split(',')))  # Create gene hypelinks
            row[-2] = ' '.join(map(snp_href, row[-2].split(';')))  # Create SNP hypelinks
            row[-1] = ' '.join(format_aa_change(row[-1].split(',')))  # Create amino acid column

            # One line can contain two or more ALT values (and consequently two or more AF values)
            # We need to split such "multiple" entries to one alt value per row
            alt_cell, afq_cell = row[3], row[4]
            for alt, afq in zip(alt_cell.split(','), afq_cell.split(',')):
                afq = float(afq)
                if afq >= af_threshold:
                    afq = '{:.3f}'.format(afq)
                    var_table.append(row[:3] + [alt] + [afq] + row[5:])
    return var_table


if __name__ == '__main__':
    args = parser.parse_args()

    # Open template and fill it with data:
    with open(args.template, 'r') as template_in, open('report.tex', 'wt') as template_out:
        stats = {}  # Container for all-reound report statistics.
        parse_covmetrics(args.covmetrics, stats)
        parse_target_pcr_metrics(args.metrics, stats)

        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)
        template = template.replace('{#SAMPLE_NAME#}', escape_latex(args.sample))
        template = template.replace('{#PANEL#}', escape_latex(args.panel))
        template = template.replace('{#AF_THRESHOLD#}', format_float(args.afthreshold, to_percent=True))
        template = template.replace('{#TOTAL_READS#}', stats['TOTAL_READS'])
        template = template.replace('{#ALIGNED_READS#}', stats['PF_UQ_READS_ALIGNED'])
        template = template.replace(
            '{#PCT_ALIGNED_READS#}', format_float(stats['PCT_PF_UQ_READS_ALIGNED'], to_percent=True))
        template = template.replace(
            '{#ALIGN_BASES_ONTARGET#}', format_float(stats['PCT_AMPLIFIED_BASES'], to_percent=True))
        template = template.replace('{#COV_MEAN#}', format_float(stats['mean_coverage']))
        template = template.replace('{#COV_20#}', format_float(stats['mean_coverage_20']))
        template = template.replace('{#COV_UNI#}', format_float(stats['coverage_uniformity'], to_percent=True))

        # Parse *.cov file: a tabular file where each row represents an amplicon.
        # Take only 5th (amplicon name) and 9th (ratio of amplicon covered) column.
        # Produce ``uncovered_amplicons`` table for amplicons with < 100% coverage (ratio below 1).
        uncovered_amplicons = []
        amplicon_count = 0
        with open(args.cov, 'r') as tfile:
            for row in csv.reader(tfile, delimiter='\t'):
                amplicon_count += 1
                cov_ratio = float(row[8])
                if cov_ratio < 1:
                    uncovered_amplicons.append([
                        escape_latex(row[4]),  # Amplicon name
                        '{:.1f}'.format(cov_ratio * 100),  # Amplicon coverage (percentage)
                    ])
        stats['ALL_AMPLICONS'] = amplicon_count
        stats['COVERED_AMPLICONS'] = amplicon_count - len(uncovered_amplicons)
        stats['RATIO_COVERED_AMPLICONS'] = (amplicon_count - len(uncovered_amplicons)) / amplicon_count
        template = template.replace('{#ALL_AMPLICONS#}', str(stats['ALL_AMPLICONS']))
        template = template.replace('{#COVERED_AMPLICONS#}', str(stats['COVERED_AMPLICONS']))
        template = template.replace('{#PCT_COV#}', format_float(stats['RATIO_COVERED_AMPLICONS'], to_percent=True))
        sentence_text = ''
        if not uncovered_amplicons:
            uncovered_amplicons = [['/', '/']]
            message = 'All amplicons are completely covered with short reads.'
            sentence_text = '\\medskip \\noindent \n {{\\boldfont ' + message + '}} \n\n \\lightfont \n'
        caption = r'Amplicons with coverage $<$ 100\%'
        table_text = list_to_tex_table(uncovered_amplicons, header=['Amplicon name', r'\% covered'], caption=caption)
        template = template.replace('{#BAD_AMPLICON_TABLE#}', sentence_text + table_text)

        # Make tables with variants for each variant caller.
        table_text = ''
        af_threshold = float(args.afthreshold)
        header_glossary = {'GEN[0].AD': 'AD', 'EFF[*].GENE': 'GENE', 'EFF[*].AA': 'AA'}
        for variant_file in args.annot_vars:
            # We have multiple (=2) variant files: this are NOT *.vcf files, but tabular files derived from them.
            # The problem is that their columns (and header names) differ depending on the variant caller tool.
            COMMON_FIELDS = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP']
            if 'gatkhc.finalvars' in variant_file.lower():
                header = COMMON_FIELDS + ['GEN[0].AD', 'FS', 'EFF[*].GENE', 'ID', 'EFF[*].AA']
                caption = 'GATK HaplotypeCaller variant calls'
            elif 'lf.finalvars' in variant_file.lower():
                header = COMMON_FIELDS + ['DP4', 'SB', 'EFF[*].GENE', 'ID', 'EFF[*].AA']
                caption = 'Lofreq variant calls'
            else:
                raise ValueError('This variant caller is not supported for report generation')

            var_table = make_variant_table(variant_file, header, af_threshold)
            header = [escape_latex(header_glossary.get(col_name, col_name)) for col_name in header]
            long_cls = [2, 3, -3, -2, -1]  # Potentially long columns are REF, ALT, GENE, ID and AA
            table_text += list_to_tex_table(
                var_table, header=header, caption=caption, long_columns=long_cls, uncut_columns=[-1])
            table_text += '\n\\newpage\n'

        template = template.replace('{#VCF_TABLES#}', table_text)

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'report.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    with open('stats.txt', 'wt') as handle:
        for key, value in sorted(stats.items()):
            print('{}\t{}'.format(key, value), file=handle)

    print("Done.")

#!/usr/bin/env python3
"""Generate amplicon multireport."""
import argparse
import math
import re
import subprocess

import numpy as np
import pandas as pd
from bokeh import palettes
from bokeh.models import HoverTool, ColumnDataSource, ColorBar, BasicTicker, LinearColorMapper, PrintfTickFormatter
from bokeh.plotting import figure, save, output_file

DECIMALS = 2

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.", nargs='+')
parser.add_argument('--covmetrics', help="Coverge metrics", nargs='+')
parser.add_argument('--cov', help="Amplicon coverage", nargs='+')
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.", nargs='+')
parser.add_argument('--vcfgatkhc', help="File with VCF GATK HaplotypeCaller data.", nargs='+')
parser.add_argument('--vcflf', help="File with VCF Lofreq data.", nargs='+')
parser.add_argument('--meancov', help="Mean amplicon coverage.", nargs='+')
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")
parser.add_argument('--afthreshold', help="Allele Frequency lower threshold.")


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


def _tsv_to_list(table_file, has_header=False, delimiter='\t', pick_columns=None):
    """Parse *.tsv file into list."""
    table_data = []
    header = None
    with open(table_file, 'r') as tfile:
        if has_header:
            header = next(tfile).strip().split(delimiter)
            if pick_columns:
                common_columns = [x for x in pick_columns if x in header]
                # Find indexes of selected columns
                temp_header = {col: i for i, col in enumerate(header)}
                pick_indexes = [temp_header[col] for col in common_columns]
                header = common_columns
        for line in tfile:
            line_content = line.strip().split(delimiter)
            if pick_columns:
                # Select only specific columns and order them as in pick_columns:
                temp_line_content = {i: col for i, col in enumerate(line_content)}
                line_content = [temp_line_content.get(i, '') for i in pick_indexes]
            table_data.append(line_content)
    return table_data, header


def cut_to_pieces(string, piece_size):
    """Cut long string into smaller pieces."""
    assert isinstance(string, str)
    string_size = len(string)
    pieces = [string[i: i + piece_size] for i in range(0, string_size, piece_size)]
    return ' '.join(pieces)


def list_to_tex_table(data, header=None, caption=None, long_columns=False, wide_columns=False, uncut_columns=[]):
    """Make a TeX table from python list."""
    lines = []
    column_alingnment = ['l' for h in header]
    if long_columns:
        for col_index in long_columns:
            column_alingnment[col_index] = 'L'
    if wide_columns:
        for col_index in wide_columns:
            column_alingnment[col_index] = 'W'

    if caption:
        lines.append('\\captionof{{table}}{{{}}}'.format(caption))

    lines.append('\\noindent')
    lines.append('\\begin{{longtable}}[l] {{ {} }}'.format('  '.join(column_alingnment)))

    if header:
        lines.append('\\rowcolor{darkblue1}')
        lines.append('\\leavevmode\\color{white}\\textbf{' +
                     '}& \\leavevmode\\color{white}\\textbf{'.join(header) + '} \\\\')
        lines.append('\\endhead')

    for line in data:
        if long_columns:
            for col_index in long_columns:
                # Insert spaces, so that wrapping can happen, for columns
                # without hyperlinks or other specified uncut columns.
                if ('\href' in line[col_index] or col_index in uncut_columns):
                    new_val = line[col_index]
                else:
                    new_val = cut_to_pieces(line[col_index], 8)
                line[col_index] = '\\multicolumn{{1}}{{m{{2.3cm}}}}{{{}}}'.format(new_val)

        if wide_columns:
            for col_index in wide_columns:
                line[col_index] = '\\multicolumn{{1}}{{m{{18cm}}}}{{{}}}'.format(line[col_index])

        lines.append(' & '.join(line) + ' \\\\')

    lines.append('\\end{longtable}')
    lines.append('{\n\\addtocounter{table}{-1}}')
    return '\n'.join(lines)


def snp_href(snpid):
    """Create LaTeX hyperlink for given SNP ID."""
    if snpid.startswith('rs'):
        url = 'http://www.ncbi.nlm.nih.gov/snp/?term={}'.format(snpid)
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


def parse_covmetrics(covmetrics_report):
    """Get coverage stats."""
    with open(covmetrics_report) as covmetrics:
        cov_data = covmetrics.readline().strip().split('\t')
        mean_coverage = float(cov_data[0])
        mean20 = float(cov_data[1])
        cov_unif = float(cov_data[2])
        return mean_coverage, mean20, cov_unif


def parse_cov(cov_file_name):
    """Parse *.cov file."""
    cov_list, _ = _tsv_to_list(cov_file_name)
    covered_amplicons = len([1 for line in cov_list if float(line[8]) >= 1])
    return cov_list, covered_amplicons


def parse_mean_covd(covd_file_name):
    """Parse *.covd file."""
    covd_list, _ = _tsv_to_list(covd_file_name)
    mean_amp_cov = {}
    for line in covd_list:
        mean_amp_cov[line[0]] = line[1]
    return mean_amp_cov


def samplelist_to_string(samplelist):
    """Convert list of lists to a string, suitable for shared variants tables."""
    samples = ['{} ({})'.format(item[0], round(float(item[1]), 3)) for item in samplelist]
    return ', '.join(samples)


def produce_warning(coverage, mean_20):
    """Produce a string with warning if coverage is less than 20% of mean."""
    return 'YES' if coverage < mean_20 else ' '


def make_heatmap(samples, variant_dict, fig_name):
    """Create a heatmap of Samples x Variants."""
    # Prepare data: keep only variants shared by at least one sample:
    shared_variants = {}
    for var_names, var_samples in variant_dict.items():
        if len(var_samples) > 1:
            shared_variants[var_names] = var_samples

    # Prepare data: make sample and variant list and NumPy array of AF.
    x_names = sorted(shared_variants.keys())
    y_names = sorted(samples)
    data = np.zeros((len(x_names), len(y_names)))
    for k, v in shared_variants.items():
        for item in v:
            variant_index = x_names.index(k)
            sample_index = y_names.index(item[0])
            data[variant_index, sample_index] = '{0:g}'.format(round(float(item[1]) * 100, DECIMALS))
    width = min(len(x_names) * 100, 1100)
    height = max(min(len(y_names) * 100, 580), 300)

    # Create DataFrame and ColumnDataSource
    df1 = pd.DataFrame(data=data, index=x_names, columns=y_names)
    df1.index.name = 'Variant'
    df1.columns.name = 'Sample'

    df = pd.DataFrame(df1.stack(), columns=['af']).reset_index()
    source = ColumnDataSource(df)

    p = figure(
        title="Shared {} across samples".format(fig_name),
        x_axis_location="above",
        tools="hover,save,pan,box_zoom,wheel_zoom,reset",
        plot_width=width,
        plot_height=height,
        x_range=x_names,
        y_range=y_names,
        toolbar_location='below',
    )
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = None
    p.xaxis.major_label_orientation = np.pi / 2

    palette = list(reversed(palettes.YlGnBu[9]))
    mapper = LinearColorMapper(palette=palette, low=np.amin(data), high=np.amax(data))
    p.rect('Variant', 'Sample', 1, 1, source=source,
           fill_color={'field': 'af', 'transform': mapper}, line_color=None, hover_line_color='black')

    p.select_one(HoverTool).tooltips = [
        ('Sample', '@Sample'),
        ('Variant', '@Variant'),
        ('Allele Frequency (%)', '@af'),
    ]

    color_bar = ColorBar(
        color_mapper=mapper,
        ticker=BasicTicker(desired_num_ticks=len(palette) + 1),
        formatter=PrintfTickFormatter(format="%d%%"),
        major_tick_line_color='black',
        major_tick_out=5, major_tick_in=0,
        label_standoff=10, border_line_color=None, location=(0, 0),
        title='Allele Frequency', title_text_font_size='7pt',
        title_standoff=5
    )

    p.add_layout(color_bar, 'right')

    output_file("{}.html".format(fig_name.replace(" ", "")), title=fig_name)
    save(p)


def format_aa_change(aa_list):
    """Create Amino Acid Change information."""
    if aa_list:
        aa = aa_list[0]
        match_obj = re.match(r'p\.([A-Za-z]*)[0-9]*([A-Za-z]*)', aa)
        if match_obj and match_obj.group(1) == match_obj.group(2):
            return 'Synon'
        else:
            return aa
    else:
        return ''


def format_float(value, decimals=DECIMALS, to_percent=False):
    """Round string representation of float to ``decimals`` decimal places.

    If ``to_percent``==True, also multiply by 100 and add latex percent sign.
    """
    value = float(value) * 100 if to_percent else float(value)
    out = '{0:g}'.format(round(value, decimals))
    if to_percent:
        out += '\\%'
    return out


def get_legacy_aa_data(varfile):
    """Get amino-acid change data from old-formatted snpeff files."""
    table, header = _tsv_to_list(varfile, has_header=True)
    aa_changes = []
    header_len = (len(header))
    for row in table:
        # Amino acid changes are written in columns following the last column with header:
        aa_changes.append(','.join(row[header_len:]))
    return aa_changes


def prepare_vcf_table(varfile, sample, variants):
    """Prepare vcf table."""
    header = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'DP', 'DP4', 'GEN[0].AD', 'SB', 'FS', 'EFF[*].GENE', 'ID',
              'EFF[*].AA']
    header_glossary = {'GEN[0].AD': 'AD', 'EFF[*].GENE': 'GENE', 'EFF[*].AA': 'AA'}

    # Parse vcf file
    vcf_table_tmp, common_columns = _tsv_to_list(varfile, has_header=True, pick_columns=header)
    if 'EFF[*].AA' not in common_columns:
        # This is for legacy data only - before snpeff process produced 'EFF[*].AA' column.
        common_columns.append('EFF[*].AA')
        legacy_aa_data = get_legacy_aa_data(varfile)
        vcf_table_tmp = [line + [legacy_aa_data[i]] for i, line in enumerate(vcf_table_tmp)]
    # Escape user inputs and change header:
    common_columns = [header_glossary[x] if (x in header_glossary) else x for x in common_columns]
    vcf_table = []
    for line_tmp in vcf_table_tmp:
        alt_cell, af_cell = line_tmp[3], line_tmp[4]
        # One line can contain two or more ALT values (and consequently two or more AF values)
        for alt, af_ in zip(alt_cell.split(','), af_cell.split(',')):
            if float(af_) >= float(args.afthreshold):
                vcf_table.append(line_tmp[:3] + [alt] + [af_] + line_tmp[5:])

    # Fill variants dict with the variants that are shared
    for line in vcf_table:
        variant_name = '{}_chr{}_{}'.format(line[-3], line[0], line[1])
        variants.setdefault(variant_name, []).append([sample, line[4]])

    # Create hyperlink, and format amino acid changes:
    vcf_table = [
        line[:-3] +
        [gene_href(line[-3])] +
        [' '.join(map(snp_href, line[-2].split(';')))] +
        [format_aa_change(line[-1].split(','))]
        for line in vcf_table
    ]
    return vcf_table, common_columns


if __name__ == '__main__':
    args = parser.parse_args()

    # Make a index list of sorted sample names (alphabetical order)
    indexlist = [i[0] for i in sorted(enumerate(args.sample), key=lambda x: x[1])]

    # Open template and fill it with data:
    with open(args.template, 'r') as template_in, open('multireport.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)
        template = template.replace('{#AF_THRESHOLD#}', format_float(args.afthreshold, to_percent=True))

        # create dict with metrics & coverage data
        d = {}
        for i, sample in enumerate(args.sample):
            d[sample] = parse_target_pcr_metrics(args.metrics[i])
            d[sample]['mean_coverage'], d[sample]['mean20'], d[sample]['cov_unif'] = \
                parse_covmetrics(args.covmetrics[i])
            # read .cov file
            d[sample]['cov_list'], d[sample]['covered_amplicons'] = parse_cov(args.cov[i])

        # make QC information table
        qc_header = ['Sample name', 'Total \\linebreak reads', 'Aligned \\linebreak reads',
                     'Aligned \\linebreak bases\\footnote{on target}', 'Mean \\linebreak coverage',
                     'Threshold \\linebreak coverage\\footnote{20\\% of mean}',
                     'Coverage \\linebreak uniformity\\footnote{\\% bases covered 20\\% above mean}',
                     'No. of \\linebreak amplicons', 'Amplicons with 100\\% coverage']
        qc_data = []
        for i in indexlist:
            sample = args.sample[i]
            qc_data.append(
                [
                    _escape_latex(sample),
                    d[sample]['TOTAL_READS'],
                    format_float(d[sample]['PCT_PF_UQ_READS_ALIGNED'], to_percent=True),
                    format_float(d[sample]['PCT_AMPLIFIED_BASES'], to_percent=True),
                    format_float(d[sample]['mean_coverage']),
                    format_float(d[sample]['mean20']),
                    format_float(d[sample]['cov_unif']) + '\\%',
                    str(len(d[sample]['cov_list'])),
                    str(d[sample]['covered_amplicons'])
                ]
            )
        qc_info = list_to_tex_table(
            qc_data, header=qc_header, long_columns=[1, 2, 3, 4, 5, 6, 7, 8], caption='QC information')
        template = template.replace('{#QCTABLE#}', qc_info)

        # Make Amplicons with coverage < 100% table
        non_cov_header = ['Sample', 'Amplicon', '\% Covered', 'Less than 20\% of mean']
        non_cov_data = []
        for i in indexlist:
            amp_cov = parse_mean_covd(args.meancov[i])
            for line in d[args.sample[i]]['cov_list']:
                if float(line[8]) < 1:
                    non_cov_data.append((
                        _escape_latex(args.sample[i]),
                        _escape_latex(line[4]),
                        '{:.1f}'.format(float(line[8]) * 100),
                        produce_warning(float(amp_cov[line[4]]), d[args.sample[i]]['mean20'])
                    ))
        if not non_cov_data:
            non_cov_data = [['/', '/']]
        non_cov_caption = 'Amplicons with coverage < 100\\%'
        non_cov_text = list_to_tex_table(non_cov_data, header=non_cov_header, caption=non_cov_caption)
        template = template.replace('{#BAD_AMPLICON_TABLE#}', non_cov_text)

        # Parse VCF files to get variant tables and dictionaries in the form {variant: samples}
        table_text = ''
        gatkhc_variants = {}
        lf_variants = {}

        for i in indexlist:

            # Make GATK variants table:
            table_text += '\\renewcommand{\\thetable}{\\arabic{table}a}'  # Set different table counting (Na)
            caption = _escape_latex('GATK HaplotypeCaller variant calls, sample {}'.format(args.sample[i]))
            vcf_table, header = prepare_vcf_table(args.vcfgatkhc[i], args.sample[i], gatkhc_variants)
            table_text += list_to_tex_table(
                vcf_table, header=header, caption=caption, long_columns=[2, 3, -2, -1], uncut_columns=[-1])
            table_text += '{\n\\addtocounter{table}{-1}}'
            table_text += '\n\\newpage\n'

            # Make LF variants table:
            table_text += '\\renewcommand{\\thetable}{\\arabic{table}b}'  # Set different table counting (Nb)
            caption = _escape_latex('Lowfreq variant calls, sample {}'.format(args.sample[i]))
            vcf_table, header = prepare_vcf_table(args.vcflf[i], args.sample[i], lf_variants)
            table_text += list_to_tex_table(
                vcf_table, header=header, caption=caption, long_columns=[2, 3, -2, -1], uncut_columns=[-1])
            table_text += '\n\\newpage\n'

        # Set table counter back to normal (N) for further tables
        table_text += '\\renewcommand{\\thetable}{\\arabic{table}}'
        template = template.replace('{#VCF_TABLES#}', table_text)

        # Crate heatmaps
        make_heatmap(args.sample, gatkhc_variants, 'GATKHC variants')
        make_heatmap(args.sample, lf_variants, 'Lowfreq variants')

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF. Call teh command two times since the first time
    args = ['pdflatex', '-interaction=nonstopmode', 'multireport.tex']
    subprocess.call(args)
    subprocess.check_call(args)

    print("Done.")

#!/usr/bin/env python2
"""Generate amplicon report."""
import argparse
import subprocess

import matplotlib as mpl
from collections import defaultdict
from matplotlib import pyplot as plt


DECIMALS = 3  # Decimal precision when presenting results:

# Set defult colors:
GENIALIS_BLUE = '#00bcd4'
GENIALIS_GREY = '#9E9E9E'
mpl.rcParams['axes.labelcolor'] = GENIALIS_GREY
mpl.rcParams['axes.edgecolor'] = GENIALIS_GREY
mpl.rcParams['lines.color'] = GENIALIS_GREY
mpl.rcParams['lines.antialiased'] = False
mpl.rcParams['text.color'] = GENIALIS_GREY
mpl.rcParams['xtick.color'] = GENIALIS_GREY
mpl.rcParams['ytick.color'] = GENIALIS_GREY

parser = argparse.ArgumentParser(description="Fill data into tex template file.")
parser.add_argument('--sample', help="Sample name.")
parser.add_argument('--panel', help="Panel name")
parser.add_argument('--cov', help="coverageBed report file (*.cov extension).")
parser.add_argument('--covd', help="coverageBed report file (*.covd extension).")
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.")
parser.add_argument('--vcf', help="File with VCF data.")
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")


def _format_latex(string):
    """Format normal string to be LaTeX compliant."""
    string = string.replace('_', '\_')
    string = string.replace('%', '\%')
    return string


def _tsv_to_list(table_file, has_header=False, delimiter='\t', pick_columns=None):
    """Parse *.tsv file into list."""
    table_data = []
    header = None
    with open(table_file, 'r') as tfile:
        if has_header:
            header = next(tfile).strip().split(delimiter)
            if pick_columns:
                # Select only specific columns and order them as in pick_columns:
                temp_header = {i: col for i, col in enumerate(header)}
                header = [temp_header[i] for i in pick_columns]
        for line in tfile:
            line_content = line.strip().split(delimiter)
            if pick_columns:
                # Select only specific columns and order them as in pick_columns:
                temp_line_content = {i: col for i, col in enumerate(line_content)}
                line_content = [temp_line_content[i] for i in pick_columns]
            table_data.append(line_content)
    return table_data, header


def parse_targetPCRmetrics(metrics_report):
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


def meanCoverage(covd):
    """Calculate Mean coverage from coverageBed (.covd) report file."""
    with open(covd) as pbc:
        cov = 0
        lines = 0
        for line in pbc:
            cov = cov + int(line.strip().split('\t')[6])
            lines = lines + 1

    return float(cov)/lines


def perAmpliconMeanCoverge(covd):
    """Calculate per amplicon mean coverage from coverageBed (.covd) report file."""
    with open(covd) as pbc:
        n_bases = defaultdict(list)
        for line in pbc:
            data = line.strip().split('\t')
            n_bases[data[4]].append(int(data[6]))

    amp_vals = {}
    for amp in n_bases:
        amp_vals[amp] = float(sum(n_bases[amp])) / len(n_bases[amp])

    return amp_vals


def coverage_uniformity(covd, mean_coverage, cutoff):
    """Parse coverageBed (.covd) report file."""
    with open(covd) as pbc:
        covered_bases = 0
        all_bases = 0
        for line in pbc:
            cov_data = line.strip().split('\t')
            if int(cov_data[6]) >= mean_coverage * cutoff:
                covered_bases = covered_bases + 1
            all_bases = all_bases + 1
    return float(covered_bases) / float(all_bases)


def make_scatterplot(x, y, hline=None, savefile=None, dpi=300):
    """Make scatterplot."""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, y, color=GENIALIS_BLUE, marker='o', linestyle=' ', markeredgewidth=0.0)
    ax.set_xlim(min(x), max(x))
    ax.set_yscale('log')
    ax.set_xlabel("Amplicon")
    ax.set_ylabel("Average coverage")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: '{:.0f}'.format(x)))
    if hline:
        ax.axhline(hline, color=GENIALIS_GREY)

    plt.savefig(savefile, dpi=dpi)


def list_to_tex_table(data, header=None, caption=None, long_columns=False):
    """Make a TeX table from python list."""
    lines = []
    column_alingnment = ['l' for h in header]
    if long_columns:
        for col_index in long_columns:
            column_alingnment[col_index] = 'L'
        lines.append('\\begin{{longtable}} {{| {} |}}'.format(' | '.join(column_alingnment)))
    else:
        lines.append('\\begin{{longtable}} {{| {} |}}'.format(' | '.join(column_alingnment)))

    if caption:
        lines.append('\\caption{{{}}}\\\\'.format(_format_latex(caption)))
    lines.append('\\hline')
    if header:
        lines.append(' & '.join(map(_format_latex, header)) + ' \\\\')
        lines.append('\\hline')
    for line in data:
        if long_columns:
            for col_index in long_columns:
                line[col_index] = '\\multicolumn{{1}}{{m{{3cm}}|}}{{{}}}'.format(line[col_index])

        lines.append(' & '.join(map(_format_latex, line)) + ' \\\\')
        lines.append('\\hline')

    lines.append('\\end{longtable}')

    return '\n'.join(lines)


if __name__ == '__main__':
    args = parser.parse_args()

    # Open template ind fill it with data:
    with open(args.template, 'r') as template_in, open('report.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)
        template = template.replace('{#SAMPLE_NAME#}', _format_latex(args.sample))
        template = template.replace('{#PANEL#}', _format_latex(args.panel))

        # Produce results:
        metrics = parse_targetPCRmetrics(args.metrics)
        pf_bases = float(metrics['PF_BASES'])
        bases_aligned = (float(metrics['ON_AMPLICON_BASES']) + float(metrics['NEAR_AMPLICON_BASES'])) / pf_bases
        mean_coverage = meanCoverage(args.covd)
        mean20 = round(0.2 * mean_coverage, DECIMALS)
        cov_unif = round(coverage_uniformity(args.covd, mean_coverage, 0.2)*100, DECIMALS)

        template = template.replace('{#TOTAL_READS#}', metrics['TOTAL_READS'])
        template = template.replace('{#ALIGNED_READS#}', metrics['PF_UQ_READS_ALIGNED'])
        template = template.replace('{#BASES_ALIGNED#}',
                                    str(round(float(metrics['PCT_AMPLIFIED_BASES'])*100, DECIMALS)))
        template = template.replace('{#BASES_TARGET#}', str(round(bases_aligned * 100, DECIMALS)))
        template = template.replace('{#COV_MEAN#}', str(round(mean_coverage, DECIMALS)))
        template = template.replace('{#COV_20#}', str(mean20))
        template = template.replace('{#COV_UNI#}', str(cov_unif))

        # Parse *.cov file:
        cov_list, _ = _tsv_to_list(args.cov)
        x = list(range(len(cov_list)))
        covered_amplicons = len([1 for line in cov_list if float(line[8]) >= 1])
        # Make a list of (amplicon_name, coverage_in_percentage) entries for uncovered amplicons:
        not_covered_amplicons = [
            (line[4], '{:.1f}'.format(float(line[8]) * 100)) for line in cov_list if float(line[8]) < 1]
        if not not_covered_amplicons:
            not_covered_amplicons = [['/', '/']]
        table_text = list_to_tex_table(not_covered_amplicons, header=['Amplicon name', '% covered'],
                                       caption='Amplicons with coverage $<$ 100 %')

        template = template.replace('{#ALL_AMPLICONS#}', str(len(cov_list)))
        template = template.replace('{#COVERED_AMPLICONS#}', str(covered_amplicons))
        template = template.replace('{#PCT_COV#}', str(round(float(covered_amplicons) / len(cov_list)*100, DECIMALS)))
        template = template.replace('{#BAD_AMPLICON_TABLE#}', table_text)

        # Scatterplot:
        perAmpCov = perAmpliconMeanCoverge(args.covd)
        cov2 = [perAmpCov[line[4]] for line in cov_list]
        make_scatterplot(x, cov2, hline=mean20, savefile='scatter.png')
        template = template.replace('{#IMAGE2#}', 'scatter.png')

        # VCF table:
        vcf_table, header = _tsv_to_list(args.vcf, has_header=True, pick_columns=[0, 1, 3, 4, 10, 2])
        # Insert space, between SNP ID's:
        vcf_table = [line[:-1] + [', '.join(line[-1].split(';'))] for line in vcf_table]
        table_text = list_to_tex_table(vcf_table, header=header, caption='ToDo!', long_columns=[5])
        template = template.replace('{#VCF_TABLE1#}', table_text)

        # Write template to 'report.tex'
        template_out.write(template)

        # link = 'http://www.ncbi.nlm.nih.gov/snp/?term=(rs1410592[Reference SNP ID]) AND homo sapiens[Organism]'.format(snp)
        # http://grch37-cancer.sanger.ac.uk/cosmic
        # http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=3602229

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'report.tex']
    subprocess.check_call(args)

    print("Done.")

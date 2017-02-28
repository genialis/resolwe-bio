#!/usr/bin/env python2
"""Generate amplicon report."""
import argparse
import subprocess

import matplotlib as mpl
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
parser.add_argument('--cov', help="coverageBed report file (*.cov extension).")
parser.add_argument('--covd', help="coverageBed report file (*.covd extension).")
parser.add_argument('--metrics', help="CollectTargetedPcrMetrics report file.")
parser.add_argument('--template', help="Report template file.")
parser.add_argument('--logo', help="Logo.")


def _format_latex(string):
    """Format normal string to be LaTeX compliant."""
    string = string.replace('_', '\_')
    string = string.replace('%', '\%')
    return string


def parse_tsv_to_list(table_file, has_header=False, delimiter='\t'):
    """Parse *.tsv file into list."""
    table_data = []
    header = None
    with open(table_file, 'r') as tfile:
        if has_header:
            header = _format_latex(next(tfile).strip()).split(delimiter)
        for line in tfile:
            line_content = _format_latex(line.strip()).split(delimiter)
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


def make_barplot(x, y, savefile=None, dpi=300):
    """Make barplot."""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.bar(x, y, 0.5, color=GENIALIS_BLUE, linewidth=0.0)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Amplicon")
    ax.set_ylabel("Per amplicon coverage")

    plt.savefig(savefile, dpi=dpi)


def make_scatterplot(x, y, hline=None, savefile=None, dpi=300):
    """Make scatterplot."""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, y, color=GENIALIS_BLUE, marker='o', linestyle=' ', markeredgewidth=0.0)
    ax.set_xlim(min(x), max(x))
    ax.set_yscale('log')
    ax.set_xlabel("Amplicon")
    ax.set_ylabel("Log total coverage")
    if hline:
        ax.axhline(hline, color=GENIALIS_GREY)

    plt.savefig(savefile, dpi=dpi)


if __name__ == '__main__':
    args = parser.parse_args()

    # Open template ind fill it with data:
    with open(args.template, 'r') as template_in, open('report.tex', 'wt') as template_out:
        template = template_in.read()
        template = template.replace('{#LOGO#}', args.logo)
        template = template.replace('{#SAMPLE_NAME#}', _format_latex(args.sample))

        # Parse metrics file:
        metrics = parse_targetPCRmetrics(args.metrics)
        template = template.replace('{#ALIGNED_READS#}', metrics['PF_UQ_READS_ALIGNED'])
        template = template.replace('{#BASES_ALIGNED#}', str(round(float(metrics['PCT_AMPLIFIED_BASES']), DECIMALS)))
        template = template.replace('{#BASES_TARGET#}', ' ')
        template = template.replace('{#COV_MEAN#}', str(round(float(metrics['MEAN_TARGET_COVERAGE']), DECIMALS)))
        mean20 = round(0.2 * float(metrics['MEAN_TARGET_COVERAGE']), DECIMALS)
        template = template.replace('{#COV_20#}', str(mean20))
        cov_unif = round(coverage_uniformity(args.covd, float(metrics['MEAN_TARGET_COVERAGE']), 0.2), DECIMALS)
        template = template.replace('{#COV_UNI#}', str(cov_unif))

        # Parse *.cov file:
        cov_list, _ = parse_tsv_to_list(args.cov)
        amplicons = [line[4] for line in cov_list]
        x = list(range(len(amplicons)))

        # Barplot:
        cov1 = [float(line[8]) for line in cov_list]
        make_barplot(x, cov1, savefile='barplot.png')
        template = template.replace('{#IMAGE1#}', 'barplot.png')

        # Scatterplot:
        cov2 = [float(line[5]) for line in cov_list]
        make_scatterplot(x, cov2, hline=mean20, savefile='scatter.png')
        template = template.replace('{#IMAGE2#}', 'scatter.png')

        # Write template to 'report.tex'
        template_out.write(template)

    # Generate PDF:
    args = ['pdflatex', '-interaction=nonstopmode', 'report.tex']
    subprocess.check_call(args)

    print("Done.")

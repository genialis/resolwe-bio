#!/usr/bin/env python3
"""Plot amplicon coverage as HTML file with Bokeh."""
import argparse
import os

import pandas as pd
from bokeh import layouts
from bokeh.embed import components
from bokeh.models import Span, HoverTool, BoxZoomTool, WheelZoomTool, ResetTool, SaveTool, PanTool, RedoTool, \
    UndoTool, Range1d
from bokeh.plotting import figure, output_file, ColumnDataSource
from jinja2 import Environment, FileSystemLoader

COLOR_CYCLE = [
    '#586e75', '#b58900', '#268bd2', '#cb4b16', '#859900', '#d33682', '#2aa198', '#dc322f', '#073642', '#6c71c4']


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot amplicon coverage as HTML file with Bokeh.')
    parser.add_argument('-i', '--infile', help='Input filename', required=True)
    parser.add_argument('-t', '--template', help='Input filename', required=True)
    parser.add_argument('-o', '--outfile', help='Output filename', required=True)
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    df = pd.read_csv(args.infile, sep=r'\s+', names=["amplicon", "meancov", "gene"])

    df['offsetcov'] = df['meancov'] + 0.1  # shift zero values by 0.1
    df = df.dropna().reset_index(drop=True)

    # Make a hover (show amplicon name on mouse-over)
    hover = HoverTool(tooltips=[("Amplicon", "@names"), ("Y value", "$y")])
    tools = [PanTool(), BoxZoomTool(), WheelZoomTool(), RedoTool(), UndoTool(), ResetTool(), hover, SaveTool()]

    # Produce plot
    output_file(args.outfile)
    fig = figure(tools=tools, width=1200, height=600, y_axis_type='log')

    # Fill plot with one point for each amplicon:
    xvals, yvals, labels, colors = [], [], [], []
    for i, (name, group) in enumerate(df.groupby('gene', sort=False)):
        xvals.extend(list(group.index))
        yvals.extend(list(group.offsetcov))
        labels.extend(list(group.amplicon))
        # Elements in the same group should have the same color. Cycle between colors in COLOR_CYCLE:
        colors.extend([COLOR_CYCLE[i % len(COLOR_CYCLE)]] * len(list(group.index)))
    data = ColumnDataSource(data=dict(x=xvals, y=yvals, names=labels, colors=colors))
    fig.circle(x='x', y='y', color='colors', size=10, source=data)

    # Make span lines on 0.05, 0.1, 0.2, 1 and 5 mutiples of mean amplicon coverage:
    mean_coverage = df.offsetcov.mean()
    span_lines = [(5.0, 'Blue'), (1.0, 'Green'), (0.2, 'Red'), (0.1, 'Purple'), (0.05, 'Magenta')]
    xmin, xmax = min(xvals) - 1, max(xvals) + 1
    for ratio, color in span_lines:
        fig.line([xmin, xmax], [mean_coverage * ratio] * 2, line_color=color, line_dash='dashed',
                 legend='{:.0f} % of mean coverage'.format(ratio * 100))

    # Customize plot:
    ymax = 2.0 * max(df.offsetcov.max() + 1000, mean_coverage * 5)
    ymin = 0.2
    fig.y_range = Range1d(ymin, ymax)
    fig.x_range = Range1d(xmin, xmax)
    fig.xaxis.major_tick_line_color = None  # Turn off x-axis major ticks
    fig.xaxis.minor_tick_line_color = None  # Turn off x-axis minor ticks
    fig.xaxis.major_label_text_font_size = '0pt'  # Hack to remove tick labels
    fig.xaxis.axis_label = 'Amplicon'
    fig.yaxis.axis_label = 'Log10 (Amplicon coverage)'

    fig.legend.location = "bottom_right"

    script, div = components(layouts.row(fig))
    with open(os.path.join(os.path.dirname(args.outfile), 'plot.js'), 'wt') as jfile:
        jfile.write('\n'.join(script.split('\n')[2:-1]))
    with open(args.outfile, 'wt') as ofile:
        env = Environment(loader=FileSystemLoader(os.path.dirname(args.template)))
        page_template = env.get_template(os.path.basename(args.template))
        html_text = page_template.render({'bokeh_div': div})  # pylint: disable=no-member
        ofile.write(html_text)


if __name__ == "__main__":
    main()

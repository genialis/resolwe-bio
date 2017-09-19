#!/usr/bin/env python3
"""Plot amplicon coverage."""
import argparse
import matplotlib.pyplot as plt
import pandas as pd

COLOR_CYCLE = [
    '#586e75', '#b58900', '#268bd2', '#cb4b16', '#859900', '#d33682', '#2aa198', '#dc322f', '#073642', '#6c71c4']


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot amplicon coverage uniformity')
    parser.add_argument('-i', '--infile', help='Input filename', required=True)
    parser.add_argument('-o', '--outfile', help='Output filename', required=True)
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    # read input file into dataframe
    df = pd.read_csv(args.infile, sep='\s+', names=["amplicon", "meancov", "gene"], )

    # shift zero values by 0.1
    df['offsetcov'] = df['meancov'] + 0.1
    df = df.dropna().reset_index(drop=True)

    # coverage stats
    mncov = df.offsetcov.mean()
    mean02 = 0.2 * mncov
    mean01 = 0.1 * mncov
    mean005 = 0.05 * mncov
    mean5 = 5.0 * mncov

    # group data by gene/ID column for coloring points
    dfgrps = df.groupby('gene', sort=False)

    xwidth = len(df.amplicon) + 2
    xrnge = range(-10, xwidth)
    ymax = 2.0 * max(df.offsetcov.max() + 1000, mean5)
    ymin = 0.2

    # Plot data
    f, ax = plt.subplots(figsize=(32, 8))  # tune figure size/aspect ratio

    for i, (name, group) in enumerate(dfgrps):
        color = COLOR_CYCLE[i % len(COLOR_CYCLE)]
        ax.plot(group.index, group.offsetcov, label=name, marker='o', ms=8, linestyle='', alpha=0.7, color=color)

    ax.plot((-2, xwidth + 2), (mean005, mean005), ls='--', c='magenta', linewidth=1.0, alpha=0.5, label="5% mean cov")
    ax.plot((-2, xwidth + 2), (mean01, mean01), ls='--', c='purple', linewidth=1.0, alpha=0.5, label="10% mean cov")
    ax.plot((-2, xwidth + 2), (mean02, mean02), 'r--', linewidth=2.0, label="20% mean cov")
    ax.plot((-2, xwidth + 2), (mean5, mean5), ls='--', c='blue', linewidth=2.0, label="5x mean cov")
    ax.plot((-2, xwidth + 2), (mncov, mncov), ls='--', c='green', linewidth=1.0, alpha=0.3)

    ax.set_xlim((-1, xwidth + 1))
    ax.set_ylim((ymin, ymax))
    ax.set_yscale('log', basey=10)
    ax.margins(0.05)
    plt.xticks([])
    plt.yticks(fontsize=18, fontweight='bold')
    plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 1.02), loc=3, mode='expand', ncol=16)

    ax.set_xlabel("amplicon", fontsize=20, fontweight='bold')
    ax.set_ylabel(r'$log_{10}$(amplicon coverage)', fontsize=20, fontweight='bold')

    # save plot to file
    f.savefig(args.outfile, bbox_inches="tight")


if __name__ == "__main__":
    main()

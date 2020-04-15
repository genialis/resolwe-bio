#!/usr/bin/python3
"""Plot the spike-ins measured abundance for samples quality control.

Scatter plot of measured abundance (in counts) of ERCC transcripts in
log2 scale and their known concentration  [attomoles/ul] in the sample
in log2 scale. The number of undetected ERCC transcripts are noted in
the legend and are marked with gray points in the scatter plot. The
coefficient of determination (R-squared value) is calculated by
comparing  measured abundance of ERCCs in the sample with  known
concentration of ERCC spike-ins.

The histogram above the scatter plot, shows the number of ERCCs per
input concentration. Missing ERCC transcripts are colored in gray.

The scatter plot helps to infer the following:

    - How well the RNA-seq experiment has measured the known
      concentration of the spike-ins. The strength of this linear
      relationship is described with coefficient of determination,
      which is the fraction of variation in the dependent variable
      (measured abundance) that is predictable from the independent
      variable (input concentration). In general, values above 0.9 for
      example are considered ideal.

    - Whether the data spans the expected dynamic range (ERCC
      abundances span a 2^20 dynamic range). Based on the limit of
      detection for the experiment, one can identify and omit any
      endogenous genes in the sample that are insufficiently expressed
      for reliable analysis.
"""

import argparse

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from resolwe_runtime_utils import warning
from scipy import stats

import utils

# File with known ERCC concetrations
ERCC_CONC_FILE = "/opt/resolwebio/assets/ERCC_table.txt"

# known SIRV concetrations
SIRV_CONC = {
    "SIRV1": 8000,
    "SIRV2": 6000,
    "SIRV3": 11000,
    "SIRV4": 7000,
    "SIRV5": 12000,
    "SIRV6": 18000,
    "SIRV7": 7000,
}

ERCC_MIXES = {
    "ercc_mix1": {"name": "concentration in Mix 1 (attomoles/ul)",},
    "ercc_mix2": {"name": "concentration in Mix 2 (attomoles/ul)",},
    "sirv_set3": {"name": "concentration in SIRV Set 3 (attomoles/ul)",},
}

# Histogram and scatter plot x-limit
XLIM = (-9, 16)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample_names", help="List of sample names.", nargs="+")
    parser.add_argument("--sample_exps", help="List of expressions.", nargs="+")
    parser.add_argument(
        "--exp_types", help="List of expression types (TPM, SMP, FPKM...)", nargs="+"
    )
    parser.add_argument(
        "--spikeins_mix", help="Used spike-ins mix.", choices=ERCC_MIXES.keys()
    )
    return parser.parse_args()


def validate_inputs(args):
    """Validate inputs."""
    # Validate that all expression types are equal.
    exp_type_set = set(args.exp_types)
    if len(exp_type_set) != 1:
        msg = "All samples should have the same expression type, but multiple expression types were given: {}."
        msg = msg.format(", ".join(exp_type_set))
        print(warning(msg))

    # Validate that same number of sample names, expression files and
    # expression types are given.
    assert len(args.sample_names) == len(args.sample_exps) == len(args.exp_types)


def get_expected(spikeins_mix, log2=False):
    """Get expected ERCC / SIRV concentrations for given spike-in mix.

    If specified, also log2 transform.
    """
    ercc_data = pd.read_csv(
        ERCC_CONC_FILE, delimiter="\t", dtype=str, index_col="ERCC ID"
    )

    if spikeins_mix in ["ercc_mix1", "ercc_mix2"]:
        series = ercc_data.loc[:, ERCC_MIXES[spikeins_mix]["name"]].astype("float")
    elif spikeins_mix == "sirv_set3":
        # This mix contains both: SIRV's and ERCC Mix 1
        sirv_series = pd.Series(SIRV_CONC)
        ercc_series = ercc_data.loc[:, ERCC_MIXES["ercc_mix1"]["name"]].astype("float")
        series = pd.concat([ercc_series, sirv_series], axis=0)
    else:
        raise ValueError("Invalid spikein mix.")

    if log2:
        series = np.log2(series)

    return series


def get_measured(
    sample_exp, sample_name, exp_type, only_zero=False, only_nonzero=False, log2=False
):
    """Get measured expression values.

    If specified, also log2 transform and only keep nonzero values.
    """
    handle = utils.gzopen(sample_exp)
    exp = pd.read_csv(handle, delimiter="\t", index_col="Gene")
    exp = exp.loc[:, "Expression"].astype("float")

    assert not (only_zero and only_nonzero)
    if only_zero:
        exp = exp[exp == 0]
    elif only_nonzero:
        exp = exp.iloc[exp.nonzero()[0]]

    if log2:
        exp = np.log2(exp)

    return exp


def merge_expected_measured(expected, measured):
    """Merge expected and measured data."""
    merged = pd.concat([expected, measured], axis=1, join_axes=[measured.index])
    return merged


def plot_histogram(histogram, all_spikeins, nonzero_spikeins, spikein_type):
    """Histogram showing the number of all and nonzero spikeins."""
    histogram.bar(all_spikeins.index, all_spikeins, color="grey", alpha=0.5, width=0.4)
    histogram.bar(
        nonzero_spikeins.index, nonzero_spikeins, color="green", alpha=0.5, width=0.4
    )

    histogram.grid(alpha=0.5, axis="y")
    histogram.set_ylabel("No. of {}s".format(spikein_type))
    histogram.set_xlim(XLIM)


def plot_scatter(scatter, zero, nonzero, exp_type):
    """Scatter plot of measured vs. expected spikein expressions."""
    scatter.scatter(nonzero.iloc[:, 0], nonzero.iloc[:, 1], alpha=0.5, color="green")

    # Nonzero values. Use linear regression to compute expected positions
    slope, intercept, r_value, _, _ = stats.linregress(
        nonzero.iloc[:, 0], nonzero.iloc[:, 1]
    )
    if not zero.empty:
        scatter.scatter(
            zero.iloc[:, 0],
            intercept + zero.iloc[:, 0] * slope,
            alpha=0.5,
            color="grey",
        )

    # Text box
    scatter.text(
        x=0.05,
        y=0.9,
        s="\n".join(
            [
                "$R^2$ = {}".format(round(r_value ** 2, 2)),
                "{} / {} transcripts not detected".format(
                    len(zero), len(zero) + len(nonzero)
                ),
            ]
        ),
        horizontalalignment="left",
        verticalalignment="top",
        transform=scatter.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5, edgecolor="gray"),
    )

    # Sytling
    scatter.grid()
    scatter.set_xlim(XLIM)
    scatter.set_xlabel("log2(Concentration [attomoles/ul])")
    scatter.set_ylabel("Sample log2({})".format(exp_type))


def plot_histogram_scatter(
    expected, zero, nonzero, spikein_type, sample_name, exp_type
):
    """Plot figure where measured vs. expected expressions are compared."""
    fig = plt.figure(figsize=(8, 6), dpi=200,)
    gspec = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 4])

    plot_histogram(
        histogram=plt.subplot(gspec[0]),
        all_spikeins=expected.value_counts(),
        nonzero_spikeins=nonzero.iloc[:, 0].value_counts(),
        spikein_type=spikein_type,
    )

    plot_scatter(
        scatter=plt.subplot(gspec[1]), zero=zero, nonzero=nonzero, exp_type=exp_type,
    )
    title = "{} ({} spike-in's)".format(sample_name, spikein_type)
    fig.suptitle(title)
    plt.savefig(title + ".png", format="png")
    plt.close()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    validate_inputs(args)

    exp_type = args.exp_types[0]
    spikeins_mix = args.spikeins_mix

    expected = get_expected(spikeins_mix, log2=True)

    min_one_has_spikeins = False  # At least one sample has spikeins = False
    warnings = []
    for sample_name, sample_exp in zip(args.sample_names, args.sample_exps):

        measured_zero = get_measured(sample_exp, sample_name, exp_type, only_zero=True)
        measured_nonzero = get_measured(
            sample_exp, sample_name, exp_type, only_nonzero=True, log2=True
        )

        merged_zero = merge_expected_measured(expected, measured_zero)
        merged_nonzero = merge_expected_measured(expected, measured_nonzero)

        # Get only ERCC spike-in's and plot the histogram-scatter figure.
        if merged_nonzero.iloc[merged_nonzero.index.str.startswith("ERCC"), :].empty:
            warnings.append(
                "All ERCC spike-ins have zero expression in sample {}".format(
                    sample_name
                )
            )
            continue

        min_one_has_spikeins = True
        plot_histogram_scatter(
            expected=expected.iloc[expected.index.str.startswith("ERCC")],
            zero=merged_zero.iloc[merged_zero.index.str.startswith("ERCC"), :],
            nonzero=merged_nonzero.iloc[merged_nonzero.index.str.startswith("ERCC"), :],
            spikein_type="ERCC",
            sample_name=sample_name,
            exp_type=exp_type,
        )

    if min_one_has_spikeins:
        for message in warnings:
            print(warning(message))
    else:
        # In case all samples have zero expression for all spikeins,
        # rather print one warning that says so (instead of printing
        # warning for each of the samples).
        print(warning("All ERCC spike-ins in all samples have zero expression."))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Expression aggregator."""

import argparse
import csv
import gzip
import json
import math

import numpy as np


def get_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Expression aggregator")
    parser.add_argument(
        "-e", "--expressions", nargs="+", help="Expressions", required=True
    )
    parser.add_argument(
        "-d", "--descriptors", nargs="+", help="Descriptors", required=True
    )
    parser.add_argument("-s", "--source", help="Source", required=True)
    parser.add_argument(
        "-t", "--expression-type", help="Expression type", required=True
    )
    parser.add_argument("-g", "--group-by", help="Group by", required=True)
    parser.add_argument("-a", "--aggregator", help="Aggregator")
    parser.add_argument("-b", "--box-plot-output", help="Box plot output file name")
    parser.add_argument(
        "-l", "--log-box-plot-output", help="Log box plot output file name"
    )
    parser.add_argument(
        "-x", "--expressions-output", help="Expressions output file name"
    )
    return parser.parse_args()


def load_expression(fn=None, sep="\t"):
    """Read expressions from file."""
    with gzip.open(fn, "rt") as f:
        reader = csv.DictReader(f, delimiter=sep)
        return {row["Gene"]: float(row["Expression"]) for row in reader}


def get_indices(descriptors, descriptor):
    """Return positions of my_element in my_list."""
    return {i for i, x in enumerate(descriptors) if x == descriptor}


def get_values(expressions, descriptors, gene, descriptor):
    """Return expressions of a gene with matching descriptor."""
    indices = get_indices(descriptors, descriptor)
    return [
        expression[gene]
        for i, expression in enumerate(expressions)
        if i in indices and gene in expression
    ]


def load_expressions(aggregator, expression_fns=[], sep="\t", descriptors=[]):
    """Read expressions from files."""
    raw_expressions = [
        load_expression(expression_fn, sep) for expression_fn in expression_fns
    ]
    if aggregator:
        raw_expressions.extend(aggregator["raw_expressions"])
        descriptors.extend(aggregator["descriptors"])
    genes = {key for raw_expression in raw_expressions for key in raw_expression.keys()}
    grouped_expressions = {
        gene: {
            descriptor: get_values(raw_expressions, descriptors, gene, descriptor)
            for descriptor in sorted(set(descriptors))
            if get_values(raw_expressions, descriptors, gene, descriptor)
        }
        for gene in genes
    }
    return raw_expressions, descriptors, grouped_expressions


def get_log_expressions(grouped_expressions):
    """Get log(expression + 1) for all expressions."""
    log_expressions = {
        gene: {
            descriptor: [
                math.log(expression + 1.0, 2.0)
                for expression in grouped_expressions[gene][descriptor]
            ]
            for descriptor in grouped_expressions[gene]
        }
        for gene in grouped_expressions
    }
    return log_expressions


def generate_statistic(expression, gene, attribute, expression_type):
    """Get box plot statistic for expressions of a single gene and attribute."""
    min_val = min(expression)
    max_val = max(expression)
    median = np.percentile(expression, 50.0)
    q1 = np.percentile(expression, 25.0)
    q3 = np.percentile(expression, 75.0)
    iqr = q3 - q1
    lowerwhisker = max(min_val, q1 - 1.5 * iqr)
    upperwhisker = min(max_val, q3 + 1.5 * iqr)
    data_count = len(expression)
    return {
        "attribute": attribute,
        "gene": gene,
        "exp_types": [expression_type],
        "min": min_val,
        "max": max_val,
        "median": median,
        "q1": q1,
        "q3": q3,
        "lowerwhisker": lowerwhisker,
        "upperwhisker": upperwhisker,
        "data_count": data_count,
    }


def get_statistics(expressions, expression_type):
    """Get box plot statistics for expressions of all genes and attributes."""
    return {
        gene: [
            generate_statistic(exp, gene, descriptor, expression_type)
            for descriptor, exp in expression.items()
        ]
        for gene, expression in expressions.items()
    }


def output_json(statistics, fname=None, compressed=False):
    """Write json to file."""
    if not compressed:
        with open(fname, "w") as f:
            json.dump(statistics, f)
    else:
        with gzip.open(fname, "wt") as f:
            json.dump(statistics, f)


def load_json(fname):
    """Read json from file."""
    with gzip.open(fname, "rt") as f:
        return json.load(f)


def check_aggregator(aggregator, source, expression_type, group_by):
    """Check aggregator fields."""
    if aggregator["source"] != source:
        raise ValueError(
            "All expressions must be annotated by the same genome database (NCBI, UCSC, ENSEMBLE,...)."
        )
    if aggregator["expression_type"] != expression_type:
        raise ValueError("All expressions must be of the same type.")
    if aggregator["group_by"] != group_by:
        raise ValueError("Group by field must be the same.")


def get_expressions_out(
    raw_expressions, descriptors, source, expression_type, group_by
):
    """Return expressions output."""
    return {
        "raw_expressions": raw_expressions,
        "descriptors": descriptors,
        "source": source,
        "expression_type": expression_type,
        "group_by": group_by,
    }


def main():
    """Compute expression statistics."""
    args = get_args()
    aggregator = None
    if args.aggregator:
        aggregator = load_json(args.aggregator)
        check_aggregator(aggregator, args.source, args.expression_type, args.group_by)
    raw_expressions, descriptors, expressions = load_expressions(
        aggregator, args.expressions, "\t", args.descriptors
    )
    log_expressions = get_log_expressions(expressions)

    if args.box_plot_output:
        statistics = get_statistics(expressions, args.expression_type)
        output_json(statistics, args.box_plot_output)
    if args.log_box_plot_output:
        log_statistics = get_statistics(log_expressions, args.expression_type)
        output_json(log_statistics, args.log_box_plot_output)
    if args.expressions_output:
        expressions_out = get_expressions_out(
            raw_expressions,
            descriptors,
            args.source,
            args.expression_type,
            args.group_by,
        )
        output_json(expressions_out, args.expressions_output, True)


if __name__ == "__main__":
    main()

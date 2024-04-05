"""Collect expression data from samples grouped by sample annotation field."""

import csv
import gzip
import json
import math

import numpy as np

from resolwe.process import (
    DataField,
    FileField,
    JsonField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


def load_expression(fn=None, sep="\t"):
    """Read expressions from file."""
    with gzip.open(fn, "rt") as f:
        reader = csv.DictReader(f, delimiter=sep)
        return {row["Gene"]: float(row["Expression"]) for row in reader}


def get_indices(descriptors, descriptor):
    """Return positions of the descriptor in the list."""
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
    open_file = gzip.open if compressed else open
    with open_file(fname, "wt") as f:
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


class ExpressionAggregator(Process):
    """Aggregate expression data based on sample metadata fields.

    The Expression aggregator process should not be run in Batch Mode, as this will create
    redundant outputs. Rather, select multiple samples below for which you wish to aggregate the
    expression matrix.
    """

    slug = "expression-aggregator"
    name = "Expression aggregator"
    process_type = "data:aggregator:expression"
    version = "1.0.0"
    category = "Quantify"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
            "storage": 100,
        },
    }
    data_name = "Expression aggregator"

    class Input:
        """Input fields to process ExpressionAggregator."""

        exps = ListField(
            DataField("expression"),
            label="Expression data",
            description="Select expression data to aggregate.",
        )
        group_by = StringField(
            label="Sample annotation field",
            description="Select sample annotation field to group by.",
        )
        expr_aggregator = DataField(
            "aggregator:expression", label="Expression aggregator", required=False
        )

    class Output:
        """Output fields to process ExpressionAggregator."""

        exp_matrix = FileField(label="Expression aggregator")
        box_plot = JsonField(label="Box plot")
        log_box_plot = JsonField(label="Log box plot")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        exp_type = StringField(label="Expression type")

    def run(self, inputs, outputs):
        """Run the analysis."""

        expression_type = inputs.exps[0].output.exp_type
        species = inputs.exps[0].output.species
        source = inputs.exps[0].output.source

        exp_matrix_output = "exp_matrix.json.gz"
        box_plot_output = "box_plot.json"
        log_box_plot_output = "log_box_plot.json"

        # sanity checks
        for exp in inputs.exps:
            if exp.output.source != source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{exp.output.source} and {source}."
                )
            if exp.output.species != species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.output.species} and {species}."
                )
            if exp.output.exp_type != inputs.exps[0].output.exp_type:
                self.error(
                    "Input samples are of different expression type: "
                    f"{exp.output.exp_type} and {expression_type}."
                )
        # check if the expression aggregator is given/valid
        aggregator = None
        if inputs.expr_aggregator:
            if inputs.expr_aggregator.output.species != inputs.exps[0].output.species:
                self.error(
                    "Cannot append expression data to the existing Expression Aggregator object that is of different species."
                    f"{inputs.expr_aggregator.output.species} is not the same as {species}."
                )

            aggregator = load_json(inputs.expr_aggregator.output.exp_matrix)

            check_aggregator(
                aggregator,
                source,
                expression_type,
                inputs.group_by,
            )

        # prepare sample annotation data
        try:
            descriptor_data = [
                e.entity.annotations[inputs.group_by] for e in inputs.exps
            ]
            if not all(descriptor_data):
                self.error(
                    f"All samples must be annotated by the selected {inputs.group_by} field."
                )
        except KeyError:
            self.error(
                f"Annotation field {inputs.group_by} is missing from one or more samples."
            )

        expression_files = [e.output.exp.path for e in inputs.exps]
        raw_expressions, descriptors, expressions = load_expressions(
            aggregator, expression_files, "\t", descriptor_data
        )
        log_expressions = get_log_expressions(expressions)

        statistics = get_statistics(expressions, expression_type)
        output_json(statistics, box_plot_output)

        log_statistics = get_statistics(log_expressions, expression_type)
        output_json(log_statistics, log_box_plot_output)

        expressions_out = get_expressions_out(
            raw_expressions,
            descriptors,
            source,
            expression_type,
            inputs.group_by,
        )
        output_json(expressions_out, exp_matrix_output, True)

        outputs.exp_matrix = exp_matrix_output
        outputs.box_plot = box_plot_output
        outputs.log_box_plot = log_box_plot_output
        outputs.exp_type = expression_type
        outputs.source = source
        outputs.species = species

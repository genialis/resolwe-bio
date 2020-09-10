"""Find genes with similar expression profile."""
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from resolwe.process import (
    DataField,
    JsonField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)


def get_expression(fname, sep="\t"):
    """Read expressions from file."""
    df = pd.read_csv(
        filepath_or_buffer=fname,
        sep=sep,
        header=0,
        index_col=0,
        compression="gzip",
        dtype={
            0: str,
            1: float,
        },
        keep_default_na=False,
    )
    df.index = df.index.map(str)
    return df


def get_mean_expression(fnames, name, sep="\t"):
    """Get mean expression for replicates of one time point."""
    dfs = [get_expression(fname, sep=sep) for fname in fnames]
    joined = pd.concat(dfs, axis=1, join="inner")
    return joined.mean(axis=1).rename(name)


def join_expressions(positions, labels, sep="\t"):
    """Join mean expressions.

    Join expressions from different time points and return only those that are
    in all samples.
    """
    dfs = []
    for position, replicates in positions:
        dfs.append(
            get_mean_expression(
                replicates,
                name=labels[position],
                sep=sep,
            )
        )

    inner = pd.concat(dfs, axis=1, join="inner")
    outer = pd.concat(dfs, axis=1, join="outer")
    excluded = sorted(outer.index.difference(inner.index))
    return inner, excluded


def calculate_spearman(x, y):
    """Calculate Spearman correlation distance between x and y."""
    return 1.0 - spearmanr(x, y).correlation


def calculate_pearson(x, y):
    """Calculate Pearson correlation distance between x and y."""
    return 1.0 - pearsonr(x, y)[0]


def is_const(values):
    """Return True, if all values are approximately equal, otherwise return False."""
    mn = np.min(values)
    mx = np.max(values)
    if mn + mx == 0.0:
        return mn == mx
    else:
        return (mx - mn) / abs(mx + mn) < 1.0e-6


def remove_const_genes(expressions):
    """Remove genes with constant expression profile across samples."""
    matches = expressions.apply(lambda x: not is_const(x), axis=1)
    return expressions.loc[matches], matches[~matches].index.tolist()


class FindSimilar(Process):
    """Find genes with similar expression profile.

    Find genes that have similar expression over time to the query gene.
    """

    slug = "find-similar"
    name = "Find similar genes"
    process_type = "data:similarexpression"
    version = "1.0.0"
    scheduling_class = SchedulingClass.INTERACTIVE
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.6.0"}},
        "resources": {"cores": 1, "memory": 8192},
        "relations": [{"type": "series"}],
    }
    data_name = "Genes similar to {{gene}}"

    class Input:
        """Input fields to process FindSimilar."""

        expressions = ListField(
            DataField("expression"),
            relation_type="series",
            label="Time series relation",
            description="Select time course to which the expressions belong to.",
        )
        gene = StringField(
            label="Query gene",
            description="Select a gene to which others are compared.",
        )
        distance = StringField(
            label="Distance metric",
            choices=[
                ("spearman", "Spearman"),
                ("pearson", "Pearson"),
            ],
            default="spearman",
        )

    class Output:
        """Output field of the process FindSimilar."""

        similar_genes = JsonField(label="Similar genes")
        source = StringField(label="Gene ID database")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def get_data_info(self, data_objects):
        """Get data time course labels and position sorted counts paths."""
        data_positions = defaultdict(list)
        labels = {}
        for data in data_objects:
            for relation in data.relations:
                if relation.type == "series":
                    position, label = next(
                        (p.position, p.label)
                        for p in relation.partitions
                        if p.entity_id == data.entity_id
                    )
                else:
                    self.error(
                        f"Relations of type series are not defined for {data.name}."
                    )

            labels[position] = label
            data_positions[position].append(data.exp.path)
        return labels, data_positions.items()

    def run(self, inputs, outputs):
        """Run the analysis."""
        for exp in inputs.expressions:
            if exp.source != inputs.expressions[0].source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{exp.source} and {inputs.expressions[0].source}."
                )
            if exp.species != inputs.expressions[0].species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.species} and {inputs.expressions[0].species}."
                )
            if exp.exp_type != inputs.expressions[0].exp_type:
                self.error(
                    "Input samples are of different Expression types: "
                    f"{exp.exp_type} and {inputs.expressions[0].exp_type}."
                )
            if exp.feature_type != inputs.expressions[0].feature_type:
                self.error(
                    "Input samples are of different Feature type: "
                    f"{exp.feature_type} and {inputs.expressions[0].feature_type}."
                )

        labels, positions = self.get_data_info(inputs.expressions)

        if len(labels) == 1:
            self.error(
                "Only one time point was provided. At least two time "
                "points are required."
            )

        expressions, excluded = join_expressions(positions, labels)
        if len(expressions.index) < 2:
            self.error("At least two genes shared across all samples are required.")

        if excluded:
            suffix = "" if len(excluded) <= 3 else ", ..."
            excluded_genes = ", ".join(excluded[:3])
            self.warning(
                "Genes not present in all of the selected samples are "
                f"excluded from the analysis. Excluded {len(excluded)} "
                f"of them ({excluded_genes + suffix})."
            )

        if inputs.gene not in expressions.index:
            self.error(
                "Selected query gene was not found. Please make sure "
                "the selected gene name can be found in all expression "
                "time courses."
            )

        expressions, removed = remove_const_genes(expressions)

        rows = expressions.index
        if len(rows) < 2:
            self.error(
                "There are less than two genes with non-constant "
                "expression across time points. Distances can not be "
                "computed."
            )

        suffix = "" if len(removed) <= 3 else ", ..."
        if removed:
            removed_genes = ", ".join(removed[:3])
            self.warning(
                f"{len(removed)} genes ({removed_genes+suffix}) have "
                "constant expression across time points. Those genes "
                "are excluded from the computation of hierarchical "
                "clustering of genes."
            )
            if inputs.gene in removed:
                self.error(
                    f"Query gene ({inputs.gene}) has constant "
                    "expression and was removed. Distances can not be "
                    "computed."
                )

        distance_map = {"pearson": calculate_pearson, "spearman": calculate_spearman}
        distance_func = distance_map[inputs.distance]
        selected_etc = expressions.loc[inputs.gene].tolist()
        expressions = expressions.drop(inputs.gene)

        distances = np.array(
            [distance_func(selected_etc, etc) for etc in expressions.as_matrix()]
        )
        genes = expressions.index
        similarity = [
            {"gene": gene, "distance": d} for gene, d in zip(genes, distances)
        ]
        similarity.sort(key=lambda x: x["distance"])
        formatted_output = {"search gene": inputs.gene, "similar genes": similarity}

        with open("similar_genes.json", "w") as f:
            json.dump(formatted_output, f)

        outputs.similar_genes = "similar_genes.json"
        outputs.source = inputs.expressions[0].source
        outputs.species = inputs.expressions[0].species
        outputs.build = inputs.expressions[0].build
        outputs.feature_type = inputs.expressions[0].feature_type

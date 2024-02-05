"""Cluster gene expression time courses."""

import json
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import spearmanr

from resolwe.process import (
    BooleanField,
    DataField,
    JsonField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def get_expression(fname, sep="\t", gene_set=[]):
    """Read expressions from file and return only expressions of genes in gene_set."""
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
    if not gene_set:
        return df
    intersection = [gene for gene in gene_set if gene in df.index]
    return df.loc[intersection]


def get_mean_expression(fnames, name, sep="\t", gene_set=[]):
    """Get mean expression for replicates of one time point."""
    dfs = [get_expression(fname, sep=sep, gene_set=gene_set) for fname in fnames]
    joined = pd.concat(dfs, axis=1, join="inner")
    return joined.mean(axis=1).rename(name)


def join_expressions(positions, labels, sep="\t", gene_set=[]):
    """Join mean expressions.

    Join expressions from different time points and return only those that are
    in all samples and the gene_set.
    """
    dfs = []
    for position, replicates in positions:
        dfs.append(
            get_mean_expression(
                replicates,
                name=labels[position],
                sep=sep,
                gene_set=gene_set,
            )
        )

    inner = pd.concat(dfs, axis=1, join="inner")
    outer = pd.concat(dfs, axis=1, join="outer")
    if gene_set:
        excluded = sorted(set(gene_set).difference(set(inner.index)))
    else:
        excluded = sorted(outer.index.difference(inner.index))
    return inner, excluded


def get_distance_metric(distance_metric):
    """Get distance metric."""
    if distance_metric == "spearman":
        return lambda x, y: 1.0 - spearmanr(x, y).correlation
    elif distance_metric == "pearson":
        return "correlation"
    return distance_metric


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


class ClusterTimeCourse(Process):
    """Cluster gene expression time courses.

    Hierarchical clustering of expression time courses.
    """

    slug = "clustering-hierarchical-etc"
    name = "Hierarchical clustering of time courses"
    process_type = "data:clustering:hierarchical:etc"
    version = "1.3.1"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {"cores": 1, "memory": 4096, "storage": 10},
        "relations": [{"type": "series"}],
    }
    data_name = "Hierarchical clustering of time courses"
    category = "Enrichment and Clustering"

    class Input:
        """Input fields to process ClusterTimeCourse."""

        expressions = ListField(
            DataField("expression"),
            relation_type="series",
            label="Time series relation",
            description="Select time course to which the expressions belong to.",
        )
        genes = ListField(
            StringField(),
            label="Gene subset",
            required=False,
            description="Select at least two genes or leave this field empty.",
        )
        gene_species = StringField(
            label="Species",
            description="Species to which the selected genes belong to. "
            "This field is required if gene subset is set.",
            required=False,
            hidden="!genes",
            allow_custom_choice=True,
            choices=[
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Homo sapiens", "Homo sapiens"),
                ("Macaca mulatta", "Macaca mulatta"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
            ],
        )
        gene_source = StringField(
            label="Gene ID database of selected genes",
            description="This field is required if gene subset is set.",
            required=False,
            hidden="!genes",
        )
        distance = StringField(
            label="Distance metric",
            choices=[
                ("euclidean", "Euclidean"),
                ("spearman", "Spearman"),
                ("pearson", "Pearson"),
            ],
            default="spearman",
        )
        linkage = StringField(
            label="Linkage method",
            choices=[
                ("single", "single"),
                ("average", "average"),
                ("complete", "complete"),
            ],
            default="average",
        )
        ordering = BooleanField(
            label="Use optimal ordering",
            description="Results in a more intuitive tree structure, "
            "but may slow down the clustering on large datasets",
            default=False,
        )

    class Output:
        """Output field of the process ClusterTimeCourse."""

        cluster = JsonField(label="Hieararhical clustering")
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
            data_positions[position].append(data.output.exp.path)
        return labels, data_positions.items()

    def get_clustering(
        self, expressions, linkage_method="average", metric="correlation", order=False
    ):
        """Compute linkage, order, and produce a dendrogram."""
        try:
            link = linkage(
                y=expressions,
                method=linkage_method,
                metric=metric,
                optimal_ordering=order,
            )
        except Exception as err:
            self.error(f"Cannot compute linkage. Original error was: {err}")
        try:
            dend = dendrogram(link, no_plot=True)
        except Exception as err:
            self.error(f"Cannot compute dendrogram. Original error was: {err}")
        return link, dend

    def run(self, inputs, outputs):
        """Run the analysis."""
        for exp in inputs.expressions:
            if exp.output.source != inputs.expressions[0].output.source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{exp.output.source} and {inputs.expressions[0].output.source}."
                )
            if exp.output.species != inputs.expressions[0].output.species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.output.species} and {inputs.expressions[0].output.species}."
                )
            if exp.output.exp_type != inputs.expressions[0].output.exp_type:
                self.error(
                    "Input samples are of different Expression types: "
                    f"{exp.output.exp_type} and {inputs.expressions[0].output.exp_type}."
                )
            if exp.output.feature_type != inputs.expressions[0].output.feature_type:
                self.error(
                    "Input samples are of different Feature type: "
                    f"{exp.output.feature_type} and {inputs.expressions[0].output.feature_type}."
                )

        if inputs.genes:
            if len(inputs.genes) == 1:
                self.error("At least two genes have to be selected.")
            if inputs.gene_species != inputs.expressions[0].output.species:
                self.error(
                    "Selected genes must belong to the same species as "
                    "expression files. Instead genes belong to "
                    f"{inputs.gene_species} while expressions belong "
                    f"to {inputs.expressions[0].output.species}."
                )
            if inputs.gene_source != inputs.expressions[0].output.source:
                self.error(
                    "Selected genes must be annotated by the same "
                    "genome database as expressions. Instead Gene IDs "
                    f"of genes are from {inputs.gene_source} and "
                    "expressions have IDs from "
                    f"{inputs.expressions[0].output.source}."
                )

        labels, positions = self.get_data_info(inputs.expressions)

        if len(labels) == 1:
            self.error(
                "Only one time point was provided. At least two time "
                "points are required to run hierarhical clustering."
            )

        expressions, excluded = join_expressions(
            positions, labels, gene_set=inputs.genes
        )
        rows = expressions.index
        if len(rows) < 2:
            if inputs.genes:
                self.error(
                    "At least two of the selected genes have to be "
                    "present in all samples to run hierarhical "
                    f"clustering. {len(rows)} found in "
                    "all samples."
                )
            else:
                self.error(
                    "At least two genes shared across all samples are "
                    "required to run hierarhical clustering. "
                    f"{len(rows)} found in all samples."
                )

        if excluded:
            suffix = "" if len(excluded) <= 3 else ", ..."
            excluded_genes = ", ".join(excluded[:3])
            self.warning(
                "Genes not present in all of the selected samples are "
                f"excluded from the analysis. Excluded {len(excluded)} "
                f"of them ({excluded_genes + suffix})."
            )

        expressions, removed = remove_const_genes(expressions)
        rows = expressions.index
        if len(rows) == 0:
            self.error(
                "All genes have constant expression across time "
                "points. Hierarchical clustering of genes cannot be "
                "computed."
            )
        if len(rows) == 1:
            self.error(
                f"Only one gene ({rows[0]}) "
                "has a non-constant expression across time points. "
                "However, hierarchical clustering of genes cannot "
                "be computed with just one gene."
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

        link, dend = self.get_clustering(
            expressions.values,
            linkage_method=inputs.linkage,
            metric=get_distance_metric(inputs.distance),
            order=inputs.ordering,
        )

        result = {
            "gene_symbols": {i: {"gene": gene} for i, gene in enumerate(rows)},
            "linkage": link.tolist(),
            "order": dend["leaves"],
        }

        with open("cluster.json", "w") as f:
            json.dump(result, f)

        outputs.cluster = "cluster.json"
        outputs.source = inputs.expressions[0].output.source
        outputs.species = inputs.expressions[0].output.species
        outputs.build = inputs.expressions[0].output.build
        outputs.feature_type = inputs.expressions[0].output.feature_type

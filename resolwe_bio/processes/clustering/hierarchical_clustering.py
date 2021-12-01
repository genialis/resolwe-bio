"""Hierarchical clustering of genes and samples."""

import json

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import spearmanr, zscore

from resolwe.process import (
    BooleanField,
    DataField,
    GroupField,
    JsonField,
    ListField,
    Persistence,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


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


def get_expressions(fnames, sep="\t", gene_set=[]):
    """Read expressions from files.

    Return only expressions of genes that are listed in all samples and in gene_set.

    """
    dfs = [get_expression(fname, sep=sep, gene_set=gene_set) for fname in fnames]
    inner = pd.concat(dfs, axis=1, join="inner")
    outer = pd.concat(dfs, axis=1, join="outer", sort=True)
    if gene_set:
        excluded = sorted(set(gene_set).difference(set(inner.index)))
    else:
        excluded = sorted(outer.index.difference(inner.index))
    return inner, excluded


def transform(expressions, log2=False, const=1.0, z_score=False, ddof=1):
    """Compute log2 and normalize expression values.

    Parameters:
    - log2: use log2(x+const) transformation
    - const: an additive constant used in computation of log2
    - z_score: use Z-score normalization
    - ddof: degrees of freedom used in computation of Z-score

    """
    if log2:
        expressions = expressions.applymap(lambda x: np.log2(x + const))
        if expressions.isnull().values.any():
            msg = "Cannot apply log2 to expression values."
            return msg
    if z_score:
        expressions = expressions.apply(
            lambda x: zscore(x, ddof=ddof), axis=1, result_type="broadcast"
        )
        expressions.fillna(value=0.0, inplace=True)
    return expressions


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


def remove_const_samples(expressions):
    """Remove samples with constant expression profile across genes."""
    matches = expressions.apply(lambda x: not is_const(x), axis=0)
    return expressions.loc[:, matches], matches.values.tolist()


def remove_const_genes(expressions):
    """Remove genes with constant expression profile across samples."""
    matches = expressions.apply(lambda x: not is_const(x), axis=1)
    return expressions.loc[matches], matches[~matches].index.tolist()


def output_json(result=dict(), fname=None):
    """Print json if fname=None else write json to file 'fname'."""
    if fname:
        with open(fname, "w") as f:
            json.dump(result, f)
    else:
        print(json.dumps({"cluster": result}, separators=(",", ":")))


class HierarchicalClusteringSamples(ProcessBio):
    """Hierarchical clustering of samples."""

    slug = "clustering-hierarchical-samples"
    name = "Hierarchical clustering of samples"
    process_type = "data:clustering:hierarchical:sample"
    version = "3.3.0"
    category = "Other"
    data_name = "Hierarchical clustering of samples"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.3.1"}
        },
        "resources": {
            "network": True,
        },
    }

    class Input:
        """Input fields to process HierarchicalClusteringSamples."""

        exps = ListField(
            DataField("expression"),
            label="Expressions",
            description="Select at least two data objects.",
        )
        advanced = BooleanField(
            label="Show advanced options",
            default=False,
        )

        class Preprocessing:
            """Preprocessing."""

            genes = ListField(
                StringField(),
                label="Gene subset",
                required=False,
                placeholder="new gene id",
                description="Select at least two genes or leave this field empty.",
            )
            source = StringField(
                label="Gene ID database of selected genes",
                description="This field is required if gene subset is set.",
                required=False,
                hidden="!preprocessing.genes",
            )
            species = StringField(
                label="Species",
                description="Species latin name. This field is required if gene subset is set.",
                allow_custom_choice=True,
                hidden="!preprocessing.genes",
                required=False,
                choices=[
                    ("Homo sapiens", "Homo sapiens"),
                    ("Mus musculus", "Mus musculus"),
                    ("Rattus norvegicus", "Rattus norvegicus"),
                    ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ],
            )
            log2 = BooleanField(
                label="Log-transform expressions",
                default=True,
                description="Transform expressions with log2(x + 1) before clustering.",
            )
            z_score = BooleanField(
                label="Z-score normalization",
                default=True,
                description="Use Z-score normalization of gene expressions before clustering.",
            )

        class Processing:
            """Processing."""

            distance_metric = StringField(
                label="Distance metric",
                default="pearson",
                choices=[
                    ("euclidean", "Euclidean"),
                    ("pearson", "Pearson"),
                    ("spearman", "Spearman"),
                ],
            )
            linkage_method = StringField(
                label="Linkage method",
                default="average",
                choices=[
                    ("single", "single"),
                    ("average", "average"),
                    ("complete", "complete"),
                ],
            )

        class Postprocessing:
            """Postprocessing."""

            order = BooleanField(
                label="Order samples optimally",
                default=True,
            )

        preprocessing = GroupField(
            Preprocessing, label="Preprocessing", hidden="!advanced"
        )
        processing = GroupField(Processing, label="Processing", hidden="!advanced")
        postprocessing = GroupField(
            Postprocessing, label="Postprocessing", hidden="!advanced"
        )

    class Output:
        """Output fields to process HierarchicalClusteringSamples."""

        cluster = JsonField(
            label="Hierarchical clustering",
            required=False,
        )

    def get_clustering_samples(
        self,
        expressions,
        distance_metric="euclidean",
        linkage_method="average",
        order=False,
    ):
        """Compute linkage, order, and produce a dendrogram."""
        try:
            link = linkage(
                y=expressions.transpose(),
                method=linkage_method,
                metric=distance_metric,
                optimal_ordering=order,
            )
        except Exception:
            self.error("Cannot compute linkage.")
        try:
            dend = dendrogram(link, no_plot=True)
        except Exception:
            self.error("Cannot compute dendrogram.")
        return link, dend

    def run(self, inputs, outputs):
        """Run analysis."""

        sample_files = []
        sample_ids = []
        sample_names = []
        gene_labels = []

        for exp in inputs.exps:

            if exp.output.source != inputs.exps[0].output.source:
                self.warning(
                    "All expression data must be annotated by the same genome database."
                )
                self.error(
                    "Sample {} has {} gene IDs, "
                    "while sample {} has {} gene IDs.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.source,
                        exp.entity.name,
                        exp.output.source,
                    )
                )

            if exp.output.species != inputs.exps[0].output.species:
                self.warning("All expressions must be of the same Species.")
                self.error(
                    "Sample {} is {}, while sample {} is {}.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.species,
                        exp.entity.name,
                        exp.output.species,
                    )
                )

            if exp.output.exp_type != inputs.exps[0].output.exp_type:
                self.warning("All expressions must be of the same Expression type.")
                self.error(
                    "Expression {} has {} expression type, "
                    "while sample {} has {} expression type.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.exp_type,
                        exp.entity.name,
                        exp.output.exp_type,
                    )
                )

            if exp.output.feature_type != inputs.exps[0].output.feature_type:
                self.warning("All expressions must be of the same Feature type.")
                self.error(
                    "Expression {} has {} feature type, "
                    "while sample {} has {} feature type.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.feature_type,
                        exp.entity.name,
                        exp.output.feature_type,
                    )
                )

            if inputs.preprocessing.genes:
                if exp.output.source != inputs.preprocessing.source:
                    self.warning(
                        "Selected genes must be annotated by the same genome database as all expression files."
                    )
                    self.error(
                        "Gene IDs are from {} database, "
                        "while sample {} has gene IDs from {} database.".format(
                            inputs.preprocessing.source,
                            exp.entity.name,
                            exp.output.source,
                        )
                    )
                if exp.output.species != inputs.preprocessing.species:
                    self.warning(
                        "Selected genes must be from the same species as all expression files."
                    )
                    self.error(
                        "Selected genes are {}, while expression {} is {}.".format(
                            inputs.preprocessing.species,
                            exp.entity.name,
                            exp.output.species,
                        )
                    )

            sample_files.append(exp.output.exp.path)
            sample_ids.append(exp.entity.id)
            sample_names.append(exp.entity.name)

        if inputs.preprocessing.genes:
            gene_labels = [gene for gene in inputs.preprocessing.genes]

        if len(gene_labels) == 1 and inputs.processing.distance_metric != "euclidean":
            self.error(
                "Select at least two genes to compute hierarchical clustering of samples with "
                "correlation distance metric or use Euclidean distance metric."
            )

        if len(sample_files) != len(sample_ids):
            self.error(
                "The number of sample files does not match the number of sample IDs."
            )

        if len(sample_files) != len(sample_names):
            self.error(
                "The number of sample files does not match the number of sample names."
            )

        if len(sample_files) < 2:
            self.error(
                "Select at least two samples to compute hierarchical clustering of samples."
            )

        expressions, excluded = get_expressions(
            fnames=sample_files, gene_set=gene_labels
        )

        if len(expressions.index) == 0:
            if not inputs.preprocessing.genes:
                self.error("The selected samples do not have any common genes.")
            else:
                self.error("None of the selected genes are present in all samples.")

        features = self.feature.filter(
            feature_id__in=list(expressions.index),
            source=inputs.exps[0].output.source,
            species=inputs.exps[0].output.species,
        )

        if (
            len(expressions.index) == 1
            and inputs.processing.distance_metric != "euclidean"
        ):
            if not inputs.preprocessing.genes:
                self.error(
                    "The selected samples contain only one common gene ({}). At least two common "
                    "genes are required to compute hierarchical clustering of samples with "
                    "correlation distance metric. Select a different set of samples or use Euclidean "
                    "distance metric.".format(
                        sorted([feature.name for feature in features])[0]
                    )
                )
            else:
                self.error(
                    "Only one of the selected genes ({}) is present in all samples but at least two "
                    "such genes are required to compute hierarchical clustering of samples with "
                    "correlation distance metric. Select more genes or use Euclidean distance "
                    "metric.".format(sorted([feature.name for feature in features])[0])
                )

        expressions = transform(
            expressions=expressions,
            log2=inputs.preprocessing.log2,
            z_score=inputs.preprocessing.z_score,
        )
        if "Cannot apply log2" in expressions:
            self.error("Cannot apply log2 to expression values.")

        if inputs.processing.distance_metric != "euclidean":
            expressions, matches = remove_const_samples(expressions)
            if len(expressions.columns) == 0:
                self.error(
                    "All of the selected samples have constant expression across genes. Hierarchical "
                    "clustering of samples cannot be computed."
                )

            if len(expressions.columns) == 1:
                samples_name = [id for i, id in enumerate(sample_names) if matches[i]][
                    0
                ]
                self.error(
                    "Only one of the selected samples ({}) has a non-constant expression across "
                    "genes. However, hierarchical clustering of samples cannot be computed with "
                    "just one sample.".format(samples_name)
                )
            removed = [name for i, name in enumerate(sample_names) if not matches[i]]
            suffix = "" if len(removed) <= 3 else ", ..."
            if removed:
                self.warning(
                    "{} of the selected samples ({}) have constant expression across genes. "
                    "Those samples are excluded from the computation of hierarchical clustering of "
                    "samples with correlation distance "
                    "metric.".format(len(removed), ", ".join(removed[:3]) + suffix)
                )

        else:
            matches = [True] * len(sample_files)

        suffix = "" if len(excluded) <= 3 else ", ..."
        if excluded:
            features = self.feature.filter(
                feature_id__in=excluded[:3],
                source=inputs.exps[0].output.source,
                species=inputs.exps[0].output.species,
            )
            excluded_names = sorted([feature.name for feature in features])

        if len(excluded) == 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    "Gene {} is present in some but not all of the selected samples. This "
                    "gene is excluded from the computation of hierarchical clustering of "
                    "samples.".format(", ".join(excluded_names))
                )
            else:
                self.warning(
                    "{} of the selected genes ({}) is missing in at least one of the selected "
                    "samples. This gene is excluded from the computation of hierarchical "
                    "clustering of samples.".format(
                        len(excluded), ", ".join(excluded_names)
                    )
                )
        if len(excluded) > 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    "{} genes ({}) are present in some but not all of the selected samples. Those "
                    "genes are excluded from the computation of hierarchical clustering of "
                    "samples.".format(len(excluded), ", ".join(excluded_names))
                )
            else:
                self.warning(
                    "{} of the selected genes ({}) are missing in at least one of the selected "
                    "samples. Those genes are excluded from the computation of hierarchical "
                    "clustering of samples.".format(
                        len(excluded), ", ".join(excluded_names)
                    )
                )

        linkage, dendrogram = self.get_clustering_samples(
            expressions,
            distance_metric=get_distance_metric(inputs.processing.distance_metric),
            linkage_method=inputs.processing.linkage_method,
            order=inputs.postprocessing.order,
        )

        sample_ids = [sample_id for i, sample_id in enumerate(sample_ids) if matches[i]]
        result = {
            "sample_ids": {
                i: {"id": sample_id} for i, sample_id in enumerate(sample_ids)
            },
            "linkage": linkage.tolist(),
            "order": dendrogram["leaves"],
        }
        output_json(result, "cluster.json")

        outputs.cluster = "cluster.json"


class HierarchicalClusteringGenes(ProcessBio):
    """Hierarchical clustering of genes."""

    slug = "clustering-hierarchical-genes"
    name = "Hierarchical clustering of genes"
    process_type = "data:clustering:hierarchical:gene"
    version = "3.3.0"
    category = "Other"
    data_name = "Hierarchical clustering of genes"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.3.1"}
        },
        "resources": {
            "network": True,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process HierarchicalClusteringGenes."""

        exps = ListField(
            DataField("expression"),
            label="Expressions",
            description="Select at least two data objects.",
        )
        advanced = BooleanField(
            label="Show advanced options",
            default=False,
        )

        class Preprocessing:
            """Preprocessing."""

            genes = ListField(
                StringField(),
                label="Gene subset",
                required=False,
                placeholder="new gene id",
                description="Select at least two genes or leave this field empty.",
            )
            source = StringField(
                label="Gene ID database of selected genes",
                description="This field is required if gene subset is set.",
                required=False,
                hidden="!preprocessing.genes",
            )
            species = StringField(
                label="Species",
                description="Species latin name. This field is required if gene subset is set.",
                allow_custom_choice=True,
                required=False,
                hidden="!preprocessing.genes",
                choices=[
                    ("Homo sapiens", "Homo sapiens"),
                    ("Mus musculus", "Mus musculus"),
                    ("Rattus norvegicus", "Rattus norvegicus"),
                    ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ],
            )
            log2 = BooleanField(
                label="Log-transform expressions",
                default=True,
                description="Transform expressions with log2(x + 1) before clustering.",
            )
            z_score = BooleanField(
                label="Z-score normalization",
                default=True,
                description="Use Z-score normalization of gene expressions before clustering.",
            )

        class Processing:
            """Processing."""

            distance_metric = StringField(
                label="Distance metric",
                default="pearson",
                choices=[
                    ("euclidean", "Euclidean"),
                    ("pearson", "Pearson"),
                    ("spearman", "Spearman"),
                ],
            )
            linkage_method = StringField(
                label="Linkage method",
                default="average",
                choices=[
                    ("single", "single"),
                    ("average", "average"),
                    ("complete", "complete"),
                ],
            )

        class Postprocessing:
            """Postprocessing."""

            order = BooleanField(
                label="Order samples optimally",
                default=True,
            )

        preprocessing = GroupField(
            Preprocessing, label="Preprocessing", hidden="!advanced"
        )
        processing = GroupField(Processing, label="Processing", hidden="!advanced")
        postprocessing = GroupField(
            Postprocessing, label="Postprocessing", hidden="!advanced"
        )

    class Output:
        """Output fields to process HierarchicalClusteringGenes."""

        cluster = JsonField(
            label="Hierarchical clustering",
            required=False,
        )

    def get_clustering_genes(
        self,
        expressions,
        distance_metric="euclidean",
        linkage_method="average",
        order=False,
    ):
        """Compute linkage, order, and produce a dendrogram."""
        try:
            link = linkage(
                y=expressions,
                method=linkage_method,
                metric=distance_metric,
                optimal_ordering=order,
            )
        except Exception:
            self.error("Cannot compute linkage.")
        try:
            dend = dendrogram(link, no_plot=True)
        except Exception:
            self.error("Cannot compute dendrogram.")
        return link, dend

    def run(self, inputs, outputs):
        """Run analysis."""

        sample_files = []
        sample_names = []
        gene_labels = []

        for exp in inputs.exps:

            if exp.output.source != inputs.exps[0].output.source:
                self.warning(
                    "All expression data must be annotated by the same genome database."
                )
                self.error(
                    "Sample {} has {} gene IDs, "
                    "while sample {} has {} gene IDs.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.source,
                        exp.entity.name,
                        exp.output.source,
                    )
                )

            if exp.output.species != inputs.exps[0].output.species:
                self.warning("All expressions must be of the same Species.")
                self.error(
                    "Sample {} is {}, while sample {} is {}.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.species,
                        exp.entity.name,
                        exp.output.species,
                    )
                )

            if exp.output.exp_type != inputs.exps[0].output.exp_type:
                self.warning("All expressions must be of the same Expression type.")
                self.error(
                    "Expression {} has {} expression type, "
                    "while sample {} has {} expression type.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.exp_type,
                        exp.entity.name,
                        exp.output.exp_type,
                    )
                )

            if exp.output.feature_type != inputs.exps[0].output.feature_type:
                self.warning("All expressions must be of the same Feature type.")
                self.error(
                    "Expression {} has {} feature type, "
                    "while sample {} has {} feature type.".format(
                        inputs.exps[0].entity.name,
                        inputs.exps[0].output.feature_type,
                        exp.entity.name,
                        exp.output.feature_type,
                    )
                )

            if inputs.preprocessing.genes:
                if exp.output.source != inputs.preprocessing.source:
                    self.warning(
                        "Selected genes must be annotated by the same genome database as all expression files."
                    )
                    self.error(
                        "Gene IDs are from {} database, "
                        "while sample {} has gene IDs from {} database.".format(
                            inputs.preprocessing.source,
                            exp.entity.name,
                            exp.output.source,
                        )
                    )
                if exp.output.species != inputs.preprocessing.species:
                    self.warning(
                        "Selected genes must be from the same species as all expression files."
                    )
                    self.error(
                        "Selected genes are {}, while expression {} is {}.".format(
                            inputs.preprocessing.species,
                            exp.entity.name,
                            exp.output.species,
                        )
                    )

            sample_files.append(exp.output.exp.path)
            sample_names.append(exp.entity.name)

        if inputs.preprocessing.genes:
            gene_labels = [gene for gene in inputs.preprocessing.genes]

        if len(sample_files) != len(sample_names):
            self.error(
                "The number of sample files does not match the number of sample names."
            )

        if len(gene_labels) == 1:
            self.error(
                "Select at least two genes to compute hierarchical clustering of genes."
            )

        if len(sample_files) == 1 and inputs.processing.distance_metric != "euclidean":
            self.error(
                "Select at least two samples to compute hierarchical clustering of genes with "
                "correlation distance metric or use Euclidean distance metric."
            )

        expressions, excluded = get_expressions(
            fnames=sample_files, gene_set=gene_labels
        )

        if len(expressions.index) == 0:
            if not inputs.preprocessing.genes:
                self.error("The selected samples do not have any common genes.")
            else:
                self.error("None of the selected genes are present in all samples.")

        features = self.feature.filter(
            feature_id__in=list(expressions.index),
            source=inputs.exps[0].output.source,
            species=inputs.exps[0].output.species,
        )

        if (
            len(expressions.index) == 1
            and inputs.processing.distance_metric != "euclidean"
        ):
            if not inputs.preprocessing.genes:
                self.error(
                    "The selected samples contain only one common gene ({}). At least two common "
                    "genes are required to compute hierarchical clustering of genes with "
                    "correlation distance metric. Select a different set of samples or use Euclidean "
                    "distance metric.".format(
                        sorted([feature.name for feature in features])[0]
                    )
                )
            else:
                self.error(
                    "Only one of the selected genes ({}) is present in all samples but at least two "
                    "such genes are required to compute hierarchical clustering of genes with "
                    "correlation distance metric. Select more genes or use Euclidean distance "
                    "metric.".format(sorted([feature.name for feature in features])[0])
                )

        expressions = transform(
            expressions=expressions,
            log2=inputs.preprocessing.log2,
            z_score=inputs.preprocessing.z_score,
        )
        if "Cannot apply log2" in expressions:
            self.error("Cannot apply log2 to expression values.")

        if inputs.processing.distance_metric != "euclidean":
            expressions, removed = remove_const_genes(expressions=expressions)
            if len(expressions.index) == 0:
                self.error(
                    "All of the selected genes have constant expression across samples. Hierarchical "
                    "clustering of genes cannot be computed."
                )

            if len(expressions.index) == 1:
                features = self.feature.filter(
                    feature_id__in=list(expressions.index),
                    source=inputs.exps[0].output.source,
                    species=inputs.exps[0].output.species,
                )
                gene_names = sorted([feature.name for feature in features])
                self.error(
                    "Only one of the selected genes ({}) has a non-constant expression across "
                    "samples. However, hierarchical clustering of genes cannot be computed with "
                    "just one gene.".format(gene_names[0])
                )
            suffix = "" if len(removed) <= 3 else ", ..."
            if removed:
                features = self.feature.filter(
                    feature_id__in=removed[:3],
                    source=inputs.exps[0].output.source,
                    species=inputs.exps[0].output.species,
                )
                removed_names = sorted([feature.name for feature in features])
                self.warning(
                    "{} of the selected genes ({}) have constant expression across samples. "
                    "Those genes are excluded from the computation of hierarchical clustering of "
                    "genes with correlation distance "
                    "metric.".format(len(removed), ", ".join(removed_names) + suffix)
                )

        suffix = "" if len(excluded) <= 3 else ", ..."
        if excluded:
            features = self.feature.filter(
                feature_id__in=excluded[:3],
                source=inputs.exps[0].output.source,
                species=inputs.exps[0].output.species,
            )
            excluded_names = sorted([feature.name for feature in features])

        if len(excluded) == 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    "Gene {} is present in some but not all of the selected samples. This "
                    "gene is excluded from the computation of hierarchical clustering of "
                    "genes.".format(", ".join(excluded_names))
                )
            else:
                self.warning(
                    "{} of the selected genes ({}) is missing in at least one of the selected "
                    "samples. This gene is excluded from the computation of hierarchical "
                    "clustering of genes.".format(
                        len(excluded), ", ".join(excluded_names)
                    )
                )
        if len(excluded) > 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    "{} genes ({}) are present in some but not all of the selected samples. Those "
                    "genes are excluded from the computation of hierarchical clustering of "
                    "genes.".format(len(excluded), ", ".join(excluded_names))
                )
            else:
                self.warning(
                    "{} of the selected genes ({}) are missing in at least one of the selected "
                    "samples. Those genes are excluded from the computation of hierarchical "
                    "clustering of genes.".format(
                        len(excluded), ", ".join(excluded_names)
                    )
                )

        linkage, dendrogram = self.get_clustering_genes(
            expressions,
            distance_metric=get_distance_metric(inputs.processing.distance_metric),
            linkage_method=inputs.processing.linkage_method,
            order=inputs.postprocessing.order,
        )

        result = {
            "gene_symbols": {
                i: {"gene": gene} for i, gene in enumerate(expressions.index)
            },
            "linkage": linkage.tolist(),
            "order": dendrogram["leaves"],
        }
        output_json(result, "cluster.json")

        outputs.cluster = "cluster.json"

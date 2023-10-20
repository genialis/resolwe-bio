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


def check_compatibility(
    exp_source,
    target_source,
    exp_species,
    target_species,
    exp_type,
    target_exp_type,
    exp_feature_type,
    target_feature_type,
    process_source,
    process_species,
    exp_name,
    target_name,
    error,
    warning,
    genes,
):
    """Check compatibility of inputs."""
    if exp_source != target_source:
        warning("All expression data must be annotated by the same genome database.")

        error(
            f"Sample {target_name} has {target_source} gene IDs, "
            f"while sample {exp_name} has {exp_source} gene IDs."
        )

    if exp_species != target_species:
        warning("All expressions must be of the same Species.")
        error(
            f"Sample {target_name} is {target_species}, while sample {exp_name} is {exp_species}."
        )

    if exp_type != target_exp_type:
        warning("All expressions must be of the same Expression type.")
        error(
            f"Expression {target_name} has {target_exp_type} expression type, "
            f"while sample {exp_name} has {exp_type} expression type."
        )

    if exp_feature_type != target_feature_type:
        warning("All expressions must be of the same Feature type.")
        error(
            f"Expression {target_name} has {target_feature_type} feature type, "
            f"while sample {exp_name} has {exp_feature_type} feature type."
        )

    if len(genes) > 0:
        if exp_source != process_source:
            warning(
                "Selected genes must be annotated by the same genome database as all expression files."
            )

            error(
                f"Gene IDs are from {process_source} database, "
                f"while sample {exp_name} has gene IDs from {exp_source} database."
            )

        if exp_species != process_species:
            warning(
                "Selected genes must be from the same species as all expression files."
            )

            error(
                f"Selected genes are {process_species}, while expression {exp_name} is {exp_species}."
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


def transform(expressions, error, log2=False, const=1.0, z_score=False, ddof=1):
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
            error("Cannot apply log2 to expression values.")

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


def get_clustering(
    expressions,
    error,
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
        error("Cannot compute linkage.")
    try:
        dend = dendrogram(link, no_plot=True)
    except Exception:
        error("Cannot compute dendrogram.")
    return link, dend


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


def output_json(result=dict(), fname="cluster.json"):
    """Print json if fname=None else write json to file 'fname'."""
    with open(fname, "w") as f:
        json.dump(result, f)


class HierarchicalClusteringSamples(ProcessBio):
    """Hierarchical clustering of samples."""

    slug = "clustering-hierarchical-samples"
    name = "Hierarchical clustering of samples"
    process_type = "data:clustering:hierarchical:sample"
    version = "3.5.2"
    category = "Enrichment and Clustering"
    data_name = "Hierarchical clustering of samples"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {"cores": 1, "memory": 4096, "storage": 10},
    }

    class Input:
        """Input fields to process HierarchicalClusteringSamples."""

        exps = ListField(
            DataField("expression"),
            label="Expressions",
            description="Select at least two data objects.",
        )

        class Preprocessing:
            """Preprocessing."""

            genes = ListField(
                StringField(),
                label="Gene subset",
                required=False,
                placeholder="new gene id, e.g. ENSG00000185982 (ENSEMBL database)",
                description="Specify at least two genes or leave this field empty.",
            )
            source = StringField(
                label="Gene ID database of selected genes",
                description="This field is required if gene subset is set, e.g. ENSEMBL, UCSC.",
                required=False,
                hidden="!preprocessing.genes",
            )
            species = StringField(
                label="Species",
                description="Specify species name. This field is required if gene subset is set.",
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

        preprocessing = GroupField(Preprocessing, label="Preprocessing")
        processing = GroupField(Processing, label="Processing")
        postprocessing = GroupField(Postprocessing, label="Postprocessing")

    class Output:
        """Output fields to process HierarchicalClusteringSamples."""

        cluster = JsonField(
            label="Hierarchical clustering",
            required=False,
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        sample_files = []
        sample_ids = []
        sample_names = []
        if inputs.preprocessing.genes:
            gene_labels = inputs.preprocessing.genes
        else:
            gene_labels = []

        for exp in inputs.exps:
            check_compatibility(
                exp_source=exp.output.source,
                target_source=inputs.exps[0].output.source,
                exp_species=exp.output.species,
                target_species=inputs.exps[0].output.species,
                exp_type=exp.output.exp_type,
                target_exp_type=inputs.exps[0].output.exp_type,
                exp_feature_type=exp.output.feature_type,
                target_feature_type=inputs.exps[0].output.feature_type,
                process_source=inputs.preprocessing.source,
                process_species=inputs.preprocessing.species,
                exp_name=exp.entity.name,
                target_name=inputs.exps[0].entity.name,
                warning=self.warning,
                error=self.error,
                genes=gene_labels,
            )

            sample_files.append(exp.output.exp.path)
            sample_ids.append(exp.entity.id)
            sample_names.append(exp.entity.name)

        if len(gene_labels) == 1 and inputs.processing.distance_metric != "euclidean":
            self.error(
                "Select at least two genes to compute hierarchical clustering of samples with "
                "correlation distance metric or use Euclidean distance metric."
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
                    "The selected samples contain only one common gene "
                    f"({[feature.name for feature in features][0]}). At least two common "
                    "genes are required to compute hierarchical clustering of samples with "
                    "correlation distance metric. Select a different set of samples or use Euclidean "
                    "distance metric."
                )
            else:
                self.error(
                    f"Only one of the selected genes ({[feature.name for feature in features][0]}) "
                    "is present in all samples but at least two such genes are required to compute "
                    "hierarchical clustering of samples with correlation distance metric. Select more "
                    "genes or use Euclidean distance metric."
                )

        expressions = transform(
            expressions=expressions,
            error=self.error,
            log2=inputs.preprocessing.log2,
            z_score=inputs.preprocessing.z_score,
            const=1.0,
        )

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
                    f"Only one of the selected samples ({samples_name}) has a non-constant expression across "
                    "genes. However, hierarchical clustering of samples cannot be computed with "
                    "just one sample."
                )
            removed = [name for i, name in enumerate(sample_names) if not matches[i]]
            if removed:
                suffix = "" if len(removed) <= 3 else ", ..."
                self.warning(
                    f"{len(removed)} of the selected samples ({', '.join(removed[:3]) + suffix}) have "
                    "constant expression across genes. Those samples are excluded from the computation "
                    "of hierarchical clustering of samples with correlation distance metric."
                )

        else:
            matches = [True] * len(sample_files)

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
                    f"Gene {excluded_names[0]} is present in some but not all of the selected samples. "
                    "This gene is excluded from the computation of hierarchical clustering of "
                    "samples."
                )
            else:
                self.warning(
                    f"{excluded} of the selected genes ({excluded_names[0]}) is missing in at least one "
                    "of the selected samples. This gene is excluded from the computation of hierarchical "
                    "clustering of samples."
                )
        if len(excluded) > 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    f"{len(excluded)} genes ({', '.join(excluded_names)}) are present in some but "
                    "not all of the selected samples. Those genes are excluded from the computation "
                    "of hierarchical clustering of samples."
                )
            else:
                self.warning(
                    f"{len(excluded)} of the selected genes ({', '.join(excluded_names)}) are missing "
                    "in at least one of the selected samples. Those genes are excluded from the "
                    "computation of hierarchical clustering of samples."
                )

        linkage, dendrogram = get_clustering(
            expressions=expressions.transpose(),
            error=self.error,
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
        output_json(result=result)

        outputs.cluster = "cluster.json"


class HierarchicalClusteringGenes(ProcessBio):
    """Hierarchical clustering of genes."""

    slug = "clustering-hierarchical-genes"
    name = "Hierarchical clustering of genes"
    process_type = "data:clustering:hierarchical:gene"
    version = "3.5.2"
    category = "Enrichment and Clustering"
    data_name = "Hierarchical clustering of genes"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {"cores": 1, "memory": 4096, "storage": 10},
    }

    class Input:
        """Input fields to process HierarchicalClusteringGenes."""

        exps = ListField(
            DataField("expression"),
            label="Expressions",
            description="Select at least two data objects.",
        )

        class Preprocessing:
            """Preprocessing."""

            genes = ListField(
                StringField(),
                label="Gene subset",
                required=False,
                placeholder="new gene id, e.g. ENSG00000185982 (ENSEMBL database)",
                description="Specify at least two genes or leave this field empty.",
            )
            source = StringField(
                label="Gene ID database of selected genes",
                description="This field is required if gene subset is set, e.g. ENSEMBL, UCSC.",
                required=False,
                hidden="!preprocessing.genes",
            )
            species = StringField(
                label="Species",
                description="Specify species name. This field is required if gene subset is set.",
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

        preprocessing = GroupField(Preprocessing, label="Preprocessing")
        processing = GroupField(Processing, label="Processing")
        postprocessing = GroupField(Postprocessing, label="Postprocessing")

    class Output:
        """Output fields to process HierarchicalClusteringGenes."""

        cluster = JsonField(
            label="Hierarchical clustering",
            required=False,
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        sample_files = []
        sample_names = []
        if inputs.preprocessing.genes:
            gene_labels = inputs.preprocessing.genes
        else:
            gene_labels = []

        for exp in inputs.exps:
            check_compatibility(
                exp_source=exp.output.source,
                target_source=inputs.exps[0].output.source,
                exp_species=exp.output.species,
                target_species=inputs.exps[0].output.species,
                exp_type=exp.output.exp_type,
                target_exp_type=inputs.exps[0].output.exp_type,
                exp_feature_type=exp.output.feature_type,
                target_feature_type=inputs.exps[0].output.feature_type,
                process_source=inputs.preprocessing.source,
                process_species=inputs.preprocessing.species,
                exp_name=exp.entity.name,
                target_name=inputs.exps[0].entity.name,
                warning=self.warning,
                error=self.error,
                genes=gene_labels,
            )

            sample_files.append(exp.output.exp.path)
            sample_names.append(exp.entity.name)

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
                    "The selected samples contain only one common gene "
                    f"({[feature.name for feature in features][0]}). At least two common "
                    "genes are required to compute hierarchical clustering of genes with "
                    "correlation distance metric. Select a different set of samples or use Euclidean "
                    "distance metric."
                )
            else:
                self.error(
                    f"Only one of the selected genes ({[feature.name for feature in features][0]}) "
                    "is present in all samples but at least two such genes are required to compute "
                    "hierarchical clustering of genes with correlation distance metric. Select more "
                    "genes or use Euclidean distance metric."
                )

        expressions = transform(
            expressions=expressions,
            error=self.error,
            log2=inputs.preprocessing.log2,
            z_score=inputs.preprocessing.z_score,
            const=1.0,
        )

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
                gene_names = [feature.name for feature in features]
                self.error(
                    f"Only one of the selected genes ({gene_names[0]}) has a non-constant expression across "
                    "samples. However, hierarchical clustering of genes cannot be computed with "
                    "just one gene."
                )

            if removed:
                suffix = "" if len(removed) <= 3 else ", ..."
                features = self.feature.filter(
                    feature_id__in=removed[:3],
                    source=inputs.exps[0].output.source,
                    species=inputs.exps[0].output.species,
                )
                removed_names = sorted([feature.name for feature in features])
                self.warning(
                    f"{len(removed)} of the selected genes ({', '.join(removed_names) + suffix}) have "
                    "constant expression across samples. Those genes are excluded from the computation "
                    "of hierarchical clustering of genes with correlation distance metric."
                )

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
                    f"Gene {excluded_names} is present in some but not all of the selected samples. "
                    "This gene is excluded from the computation of hierarchical clustering of "
                    "genes."
                )
            else:
                self.warning(
                    f"{len(excluded)} of the selected genes ({excluded_names[0]}) is missing in at least "
                    "one of the selected samples. This gene is excluded from the computation of "
                    "hierarchical clustering of genes."
                )
        if len(excluded) > 1:
            if not inputs.preprocessing.genes:
                self.warning(
                    f"{len(excluded)} genes ({', '.join(excluded_names)}) are present in some but "
                    "not all of the selected samples. Those genes are excluded from the computation "
                    "of hierarchical clustering of genes."
                )
            else:
                self.warning(
                    f"{len(excluded)} of the selected genes ({', '.join(excluded_names)}) are "
                    "missing in at least one of the selected samples. Those genes are excluded from "
                    "the computation of hierarchical clustering of genes."
                )

        linkage, dendrogram = get_clustering(
            expressions=expressions,
            error=self.error,
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
        output_json(result=result)

        outputs.cluster = "cluster.json"

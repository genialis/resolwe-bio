"""Beta principal component analysis process."""

import json

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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


def component_top_factors(component, allgenes_array, max_size=20):
    """Return top 20 absolute factors."""
    abs_component = np.abs(component)
    size = min(component.size, max_size)
    unordered_ixs = np.argpartition(abs_component, -size)[-size:]
    ixs = unordered_ixs[np.argsort(abs_component[unordered_ixs])[::-1]]
    if ixs.size == 0:
        return []
    return list(zip(np.array(allgenes_array)[ixs].tolist(), component[ixs].tolist()))


def get_pca(expressions=pd.DataFrame(), gene_labels=[]):
    """Compute PCA."""
    if not gene_labels:
        gene_labels = expressions.index
    skipped_gene_labels = list(set(gene_labels).difference(expressions.index))

    if expressions.shape[0] < 2 or expressions.shape[1] < 2:
        coordinates = [[0.0, 0.0] for _ in range(expressions.shape[1])]
        all_components = [[], []]
        all_explained_variance_ratios = [0.0, 0.0]
    else:
        pca = PCA(n_components=2, whiten=True)
        pca_expressions = pca.fit_transform(expressions.transpose())

        coordinates = [
            t[:2].tolist() if len(t) > 1 else [t[0], 0.0] for t in pca_expressions
        ]
        all_components = [
            component_top_factors(component=component, allgenes_array=gene_labels)
            for component in pca.components_
        ]
        if np.isnan(pca.explained_variance_ratio_).any():
            all_explained_variance_ratios = [0.0 for _ in pca.explained_variance_ratio_]
        else:
            all_explained_variance_ratios = pca.explained_variance_ratio_.tolist()

    result = {
        "coordinates": coordinates,
        "all_components": all_components,
        "all_explained_variance_ratios": all_explained_variance_ratios,
        "skipped_gene_labels": skipped_gene_labels,
        "warning": None,
    }

    return result


def save_pca(result={}, sample_ids=[], max_size=10):
    """Save PCA."""
    output_fn = "pca.json"
    data = {
        "flot": {
            "data": result["coordinates"],
            "xlabel": "PC 1",
            "ylabel": "PC 2",
            "sample_ids": sample_ids,
        },
        "zero_gene_symbols": result["skipped_gene_labels"],
        "components": result["all_components"][:max_size],
        "all_components": result["all_components"],
        "explained_variance_ratios": result["all_explained_variance_ratios"][:max_size],
        "all_explained_variance_ratios": result["all_explained_variance_ratios"],
    }

    with open(output_fn, "w") as outfile:
        json.dump(data, outfile, separators=(",", ":"), allow_nan=False)


def read_csv(fname):
    """Read CSV file and return Pandas DataFrame."""
    csv = pd.read_csv(
        filepath_or_buffer=fname,
        sep="\t",
        header=0,
        index_col=0,
        dtype={
            0: str,
            1: float,
        },
        keep_default_na=False,
    )
    csv.index = csv.index.map(str)
    return csv


def transform(expressions, log2=True, const=1.0, z_score=True, error=None):
    """Compute log2 and normalize expression values.

    Parameters:
    - log2: use log2(x+const) transformation
    - const: an additive constant used in computation of log2
    - z_score: use Z-score normalization using StandardScaler
    """
    if log2:
        expressions = expressions.applymap(lambda x: np.log2(x + const))
        if expressions.isnull().values.any():
            error(
                "Cannot apply log2 to expression values. Some expression values are negative."
            )
    if z_score:
        try:
            scaler = StandardScaler()
            expressions = expressions.transpose()
            expressions_arr = scaler.fit_transform(expressions)
            expressions = pd.DataFrame(
                expressions_arr, index=expressions.index, columns=expressions.columns
            ).transpose()
        except ValueError:
            error("Cannot apply Z-score normalization to expression values.")

    return expressions


def get_csv(fnames):
    """Read CSV files, perform transformations and return Pandas DataFrame."""
    expressions = [read_csv(fname) for fname in fnames]
    return pd.concat(expressions, axis=1, join="inner")


def run_pca(
    exp_paths,
    rc_paths,
    sample_ids,
    gene_labels,
    log2_transform,
    standard_scaler,
    error,
):
    """Read data, run PCA, and output results."""
    expressions = get_csv(exp_paths)
    # Filter out genes with low expression values
    if rc_paths:
        read_counts = get_csv(rc_paths)
        percentile_80 = np.percentile(read_counts, 80, axis=1)
        retained_genes = percentile_80 > 100
        expressions = expressions.loc[retained_genes]

    if expressions.empty:
        error(
            "Gene selection and filtering resulted in no genes. Please select different samples or genes."
        )

    expressions = transform(
        error=error,
        expressions=expressions,
        log2=log2_transform,
        z_score=standard_scaler,
    )

    if gene_labels:
        gene_labels = list(set(gene_labels).intersection(expressions.index))
        gene_labels.sort()
        expressions = expressions.loc[gene_labels]

    result = get_pca(expressions=expressions, gene_labels=gene_labels)
    save_pca(result, sample_ids)


class PrinicipalComponentAnalysis(Process):
    """Principal component analysis process (beta)."""

    slug = "pca-beta"
    name = "Principal component analysis (beta)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/common:4.1.1",
            },
        },
        "resources": {
            "cores": 1,
            "memory": 4096,
            "storage": 10,
        },
    }
    data_name = "{{ alignment|name|default('?') }}"
    version = "0.1.0"
    process_type = "data:pca"
    category = "Enrichment and Clustering"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.TEMP

    description = "Principal component analysis (PCA)"

    class Input:
        """Input fields."""

        exps = ListField(DataField(data_type="expression"), label="Expressions")
        genes = ListField(
            StringField(),
            label="Gene subset",
            required=False,
            placeholder="new gene id, e.g. ENSG00000185982 (ENSEMBL database)",
            description="Specify at least two genes or leave this field empty. "
            "To subset your analysis to only a subset of genes, specify them here.",
        )
        source = StringField(
            label="Gene ID database of selected genes",
            description="This field is required if gene subset is set.",
            required=False,
            hidden="!genes",
        )
        species = StringField(
            label="Species",
            description="Species latin name. This field is required if gene subset is set.",
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
            ],
            required=False,
            hidden="!genes",
        )
        low_expression_filter = BooleanField(
            label="Low expression filter",
            default=True,
            description="Filter the input expression matrix to remove genes "
            "with low expression values. Only genes with 80th percentile above "
            "read count of 100 are retained. "
            "This option only works if raw counts are available.",
        )
        log2 = BooleanField(
            label="Log-transform expressions",
            default=True,
            description="Transform expressions with log2(x + 1) before performing PCA.",
        )
        standard_scaler = BooleanField(
            label="Transform input data using StandardScaler",
            default=True,
            description="Apply the StandardScaler transformation before performing PCA.",
        )

    class Output:
        """Output fields."""

        pca = JsonField(label="PCA")

    def run(self, inputs, outputs):
        """Run the process."""

        expressions = inputs.exps

        if inputs.genes and not (inputs.source and inputs.species):
            self.error(
                "Gene ID database and Species must be provided if gene subset is selected."
            )

        for exp in expressions:
            if exp.output.source != expressions[0].output.source:
                self.error(
                    "Input samples are of different Gene ID databases: "
                    f"{exp.output.source} and {expressions[0].output.source}."
                )
            if exp.output.species != expressions[0].output.species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.output.species} and {expressions[0].output.species}."
                )
            if exp.output.build != expressions[0].output.build:
                self.error(
                    "Input samples are of different Build: "
                    f"{exp.output.build} and {expressions[0].output.build}."
                )
            if exp.output.feature_type != expressions[0].output.feature_type:
                self.error(
                    "Input samples are of different Feature type: "
                    f"{exp.output.feature_type} and {expressions[0].output.feature_type}."
                )

            if inputs.genes:
                if exp.output.source != inputs.source:
                    self.error(
                        "Selected genes must be annotated by the same genome database as all expression files."
                        f"Gene IDs are from {inputs.source} database while sample {exp.entity.name} has gene IDs from {exp.output.source} database."
                    )
                if exp.output.species != inputs.species:
                    self.error(
                        "Selected genes must be from the same species as all expression files."
                        f"Selected genes are from {inputs.species}, while expression from sample {exp.entity.name} is from {exp.output.species}."
                    )

        expression_paths = [exp.output.exp.path for exp in expressions]
        sample_ids = [exp.entity.id for exp in expressions]

        if inputs.low_expression_filter:
            try:
                rc_paths = [exp.output.rc.path for exp in expressions]
            except AttributeError:
                rc_paths = None
                self.warning(
                    "Cannot apply low expression filter to input samples. Not all samples have raw counts available."
                )
        else:
            rc_paths = None

        try:
            run_pca(
                exp_paths=expression_paths,
                rc_paths=rc_paths,
                sample_ids=sample_ids,
                gene_labels=inputs.genes,
                log2_transform=inputs.log2,
                standard_scaler=inputs.standard_scaler,
                error=self.error,
            )
        except Exception as e:
            self.error(f"PCA computation failed: {e}")

        outputs.pca = "pca.json"

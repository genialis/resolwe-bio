#!/usr/bin/env python3
"""Principal components analysis."""

import argparse
import json

import numpy as np
import pandas as pd
from resolwe_runtime_utils import warning
from sklearn.decomposition import PCA


def get_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="PCA")
    parser.add_argument(
        "--sample-files", "-f", nargs="+", help="Sample file names", required=True
    )
    parser.add_argument(
        "--sample-ids", "-i", nargs="+", help="Sample IDs", required=True
    )
    parser.add_argument("--gene-labels", "-g", nargs="+", help="Filter genes by label")
    parser.add_argument(
        "--components", "-c", help="Number of PCA components", type=int, default=2
    )
    parser.add_argument("--output-fn", "-o", help="Output file name")
    return parser.parse_args()


def component_top_factors(component, allgenes_array, max_size=20):
    """Return top 20 absolute factors."""
    abs_component = np.abs(component)
    size = min(component.size, max_size)
    unordered_ixs = np.argpartition(abs_component, -size)[-size:]
    ixs = unordered_ixs[np.argsort(abs_component[unordered_ixs])[::-1]]
    if ixs.size == 0:
        return []
    return list(zip(np.array(allgenes_array)[ixs].tolist(), component[ixs].tolist()))


def get_pca(expressions=pd.DataFrame(), n_components=2, gene_labels=[]):
    """Compute PCA."""
    if not gene_labels:
        gene_labels = expressions.index
    skipped_gene_labels = list(set(gene_labels).difference(expressions.index))

    if expressions.shape[0] < 2 or expressions.shape[1] < 2:
        coordinates = [[0.0, 0.0] for i in range(expressions.shape[1])]
        all_components = [[], []]
        all_explained_variance_ratios = [0.0, 0.0]
    else:
        pca = PCA(n_components=n_components, whiten=True)
        pca_expressions = pca.fit_transform(expressions.transpose())

        coordinates = [
            t[:2].tolist() if len(t) > 1 else [t[0], 0.0] for t in pca_expressions
        ]
        all_components = [
            component_top_factors(component, gene_labels)
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

    if expressions.empty:
        print(
            warning(
                "Gene selection and filtering resulted in no genes. Please select different samples or genes."
            )
        )

    return result


def save_pca(result={}, sample_ids=[], output_fn=None, max_size=10):
    """Save PCA."""
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

    if output_fn:
        with open(output_fn, "w") as outfile:
            json.dump(data, outfile, separators=(",", ":"), allow_nan=False)
    else:
        print(json.dumps(data, separators=(",", ":"), allow_nan=False))


def read_csv(fname):
    """Read CSV file and return Pandas DataFrame."""
    csv = pd.read_csv(
        filepath_or_buffer=fname,
        sep="\t",
        header=0,
        index_col=0,
        dtype={0: str, 1: float,},
        keep_default_na=False,
    )
    csv.index = csv.index.map(str)
    return csv


def get_csv(fnames):
    """Read CSV files and return Pandas DataFrame."""
    expressions = [read_csv(fname) for fname in fnames]
    return pd.concat(expressions, axis=1, join="inner")


def main():
    """Read data, run PCA, and output results."""
    args = get_args()
    expressions = get_csv(args.sample_files)

    if args.gene_labels:
        gene_labels = set(args.gene_labels).intersection(expressions.index)
        expressions = expressions.loc[gene_labels]

    result = get_pca(expressions, args.components, args.gene_labels)
    save_pca(result, args.sample_ids, args.output_fn)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create gene set table."""
import argparse
from pathlib import Path

import pandas as pd
from resolwe_runtime_utils import send_message, warning


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create gene sets from DGE results table."
    )
    parser.add_argument(
        "--dge_file", required=True, help="Differential expressions file."
    )
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument(
        "--analysis_name", type=str, required=True, help="Analysis name."
    )
    parser.add_argument("--tool", type=str, required=True, help="Tool name.")
    parser.add_argument("--logfc", type=float, default=1, help="LogFC threshold.")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold")
    return parser.parse_args()


def save_genes(data_frame, outfname):
    """Save gene ids."""
    data_frame["gene_id"].to_csv(
        outfname,
        header=False,
        index=False,
        sep="\n",
        compression="gzip",
    )
    return outfname


def generate_name(analysis_name, tool_name, logfc=1, fdr=0.05):
    """Generate name for gene sets."""
    if " " in analysis_name:
        analysis_name = analysis_name.replace(" ", "_")
    if " " in tool_name:
        tool_name = tool_name.replace(" ", "_")
    return f"{analysis_name}_{tool_name}_logFC{logfc}_FDR{fdr}"


def create_gene_sets(exp_file, logfc=1, fdr=0.05):
    """Create gene sets from DGE results table.

    Create a gene set file with up-regulated (log2FC > 1, fdr < 0.05),
    down-regulated (log2FC > -1, fdr < 0.05) and all genes.
    """
    df = pd.read_csv(
        filepath_or_buffer=exp_file,
        sep="\t",
        header=0,
        compression="gzip",
        keep_default_na=False,
    )

    up = df.loc[(df["logfc"] > logfc) & (df["fdr"] < fdr)]
    down = df.loc[(df["logfc"] < -logfc) & (df["fdr"] < fdr)]
    return {"up": up, "down": down, "all": df}


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()
    gene_sets = create_gene_sets(args.dge_file, args.logfc, args.fdr)

    fname_prefix = generate_name(args.analysis_name, args.tool, args.logfc, args.fdr)

    out_dir = Path(args.out_dir)
    if not out_dir.exists():
        out_dir.mkdir()

    for name, data in gene_sets.items():
        if data.empty:
            send_message(
                warning(f"No {name}-regulated genes. Gene set was not created.")
            )
        else:
            save_genes(data, out_dir / f"{fname_prefix}_{name}.tab.gz")


if __name__ == "__main__":
    main()

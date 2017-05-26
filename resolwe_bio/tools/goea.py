#!/usr/bin/env python3
"""Run gene ontology enrichment analysis."""
from subprocess import Popen, PIPE, DEVNULL

import argparse
import tempfile
import resdk


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run gene ontology enrichment analysis.")
    parser.add_argument('source_db', type=str, help="Source database.")
    parser.add_argument('target_db', type=str, help="Target database.")
    parser.add_argument('species', type=str, help="Species.")
    parser.add_argument('feature_ids', help="Feature IDs file")
    parser.add_argument('obo', help="Ontology object file")
    parser.add_argument('gaf', help="GAF annotation file")
    parser.add_argument('--pval', type=float, default=0.1, help="P-value threshold")
    parser.add_argument('--min_genes', type=int, default=1, help="Minimum number of genes on a GO term.")
    return parser.parse_args()


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    res = resdk.Resolwe()

    with open(args.feature_ids) as gene_file:
        genes = [gene.strip() for gene in gene_file]

    org_features = res.feature.filter(source=args.source_db, feature_id=genes)

    if len(org_features) == 0:
        print('{"proc.error":"No genes were fetched from the knowledge base."}')
        exit(1)

    species = set(feature.species for feature in org_features)

    if len(species) != 1:
        print('{"proc.error":"Input genes belong to multiple species."}')
        exit(1)
    else:
        species = species.pop()

    if args.species == species and args.source_db == args.target_db:
        target_ids = genes
    else:
        features = res.mapping.filter(source_db=args.source_db, target_db=args.target_db, source_id=genes)

        if len(features) == 0:
            print('{"proc.error":"Failed to map features."}')
            exit(1)

        target_ids = [str(feature.target_id) for feature in features]

        if len(genes) > len(target_ids):
            print('{"proc.warning":"Not all features could be mapped."}')

    with tempfile.NamedTemporaryFile() as input_genes:
        input_genes.write(' '.join(target_ids).encode("UTF-8"))
        input_genes.flush()
        process = Popen(['processor', str(args.pval), str(args.min_genes), args.obo, args.gaf, input_genes.name],
                        stdout=PIPE,
                        stderr=DEVNULL
                        )
        out, err = process.communicate()

        with open('terms.json', 'w') as f:
            f.write(out.decode("UTF-8"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Build gene info (mm10)."""
import argparse
from . import build_gene_info_common as build


parser = argparse.ArgumentParser(description='Build gene info (mm10).')
parser.add_argument('--annotation', help='Annotation (GTF) file')
parser.add_argument('--gene_info', help='NCBI "Gene_Info" file')
parser.add_argument('--output', help='Output "GeneInfo" file')

args = parser.parse_args()

all_gene_ids = set()
geneId2GeneName = {}
GeneID2EntrezID = {}
GeneID2synonyms = {}
GeneID2description = {}
GeneID2mgi = {}
GeneID2Ensembl = {}

with open(args.annotation) as annotation_file:
    for line in annotation_file:
        if "gene_id" in line:
            gene_id = build.get_gene_id(line)
            all_gene_ids.add(gene_id)
            geneId2GeneName[gene_id] = build.get_gene_name(line)

with open(args.gene_info) as gene_info:
    gene_info.readline()

    for line in gene_info:
        if "MGI:" in line and "Ensembl" in line:
            a = line.split('\t')
            ens_id = build.get_ensembl_id(a[5])
            mgi_id = build.get_mgi_id(a[5])

            GeneID2synonyms[a[2]] = a[4]
            GeneID2description[a[2]] = a[8]
            GeneID2EntrezID[a[2]] = a[1]
            GeneID2mgi[a[2]] = mgi_id
            GeneID2Ensembl[a[2]] = ens_id


with open(args.output, "w") as f:
    f.write('\t'.join(["Gene ID", "Gene Name", "Synonyms", "Gene Products", "Entrez ID",
                       "MGI ID", "Ensembl ID"]) + '\n')

    for gene_id in all_gene_ids:
        if gene_id not in GeneID2EntrezID:
            continue

        entrez_id = GeneID2EntrezID[gene_id]
        gene_name = geneId2GeneName[gene_id]

        synonyms = GeneID2synonyms[gene_id] if gene_id in GeneID2synonyms else 'N/A'
        description = GeneID2description[gene_id] if gene_id in GeneID2description else 'N/A'
        mgi_id = GeneID2mgi[gene_id] if gene_id in GeneID2mgi else 'N/A'
        ensembl_id = GeneID2Ensembl[gene_id] if gene_id in GeneID2Ensembl else 'N/A'

        f.write('\t'.join([gene_id, gene_name, synonyms, description, entrez_id,
                           mgi_id, ensembl_id]) + '\n')

#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Build gene info (Homo sapiens)."""
import argparse
from . import build_gene_info_common as build

parser = argparse.ArgumentParser(description='Build gene info (Homo sapiens).')
parser.add_argument('--annotation', help='Annotation (GTF) file')
parser.add_argument('--gene_info', help='NCBI "Gene_Info" file')
parser.add_argument('--uniprotKB', help='Uniprot mapping file')
parser.add_argument('--output', help='Output "GeneInfo" file')

args = parser.parse_args()

all_gene_ids = set()
geneId2GeneName = {}
GeneID2EntrezID = {}
GeneID2synonyms = {}
GeneID2description = {}
GeneID2Omim = {}
EntrezID2UniprotKB = {}

with open(args.annotation) as annotation_file:
    for line in annotation_file:
        if "gene_id" in line:
            gene_id = build.get_gene_id(line)
            all_gene_ids.add(gene_id)
            geneId2GeneName[gene_id] = build.get_gene_name(line)

with open(args.gene_info) as gene_info:
    gene_info.readline()

    for line in gene_info:
        if "Ensembl" in line:
            a = line.split('\t')
            ens_id = build.get_ensembl_id(a[5])

            GeneID2synonyms[ens_id] = a[4]
            GeneID2description[ens_id] = a[8]
            GeneID2EntrezID[ens_id] = a[1]
            GeneID2Omim[ens_id] = build.get_omim_id(a[5])

with open(args.uniprotKB) as uniprot:
    for line in uniprot:
        line = line.split('\t')
        if line[2] != '':
            EntrezID2UniprotKB[line[2]] = line[0]

with open(args.output, "w") as f:
    f.write('\t'.join(["Gene ID", "Gene Name", "Synonyms", "Gene Products", "Entrez ID",
                       "OMIM ID", "UniprotKB ID"]) + '\n')

    for gene_id in all_gene_ids:
        if gene_id not in GeneID2EntrezID:
            continue

        entrez_id = GeneID2EntrezID[gene_id]
        gene_name = geneId2GeneName[gene_id]

        synonyms = GeneID2synonyms[gene_id] if gene_id in GeneID2synonyms else 'N/A'
        description = GeneID2description[gene_id] if gene_id in GeneID2description else 'N/A'
        omim_id = GeneID2Omim[gene_id] if gene_id in GeneID2Omim else 'N/A'

        uniprotKB = 'N/A'
        if gene_id in GeneID2EntrezID and GeneID2EntrezID[gene_id] in EntrezID2UniprotKB:
            uniprotKB = EntrezID2UniprotKB[GeneID2EntrezID[gene_id]]

        f.write('\t'.join([gene_id, gene_name, synonyms, description, entrez_id,
                           omim_id, uniprotKB]) + '\n')

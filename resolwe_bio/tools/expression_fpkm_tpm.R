#!/usr/bin/Rscript
require('rtracklayer') ; require('argparse')

parser = ArgumentParser(description='Calculate expression values in FPKM and TPM units.')
parser$add_argument('annotation', help='Gene model (GTF/GFF) annotation file')
parser$add_argument('format', help='Annotation version/format')
parser$add_argument('feature_type', help='Feature type (e.g. exon)')
parser$add_argument('id_attribute', help='GFF attribute to be used as feature ID (e.g. gene_id)')
parser$add_argument('counts_file', help='Read-count data')

args = parser$parse_args(commandArgs(trailingOnly=TRUE))

# retrieve exon lengths
annot <- import.gff(args$annotation, format=args$format, genome="NA", feature.type=args$feature_type)
exon_lengths=width(annot)
names(exon_lengths)=elementMetadata(annot)[,args$id_attribute]

# calculate gene lengths
exon_lengths_by_gene=split(exon_lengths, names(exon_lengths))
gene_lengths=sapply(exon_lengths_by_gene, sum)

# normalise counts by library size
counts <- read.delim(args$counts_file, header = FALSE)
rownames(counts) <- counts[,1]
counts[,1] <- NULL
lib_size=colSums(counts)
ncounts=t(t(counts)/lib_size)

# normalise counts by gene length
common_genes=intersect(row.names(ncounts), names(gene_lengths))
subset_ncounts=ncounts[row.names(ncounts) %in% common_genes,]
gene_lengths=gene_lengths[names(gene_lengths) %in% common_genes]
ncounts=subset_ncounts/gene_lengths

# obtain RPKMs
rpkm=ncounts*1e9

# obtain TPM
tpm = exp(log(rpkm) - log(sum(rpkm)) + log(1e6))

# write data tables
write.table(rpkm, file="fpkm_wo_header.tab", sep="\t", na="0", quote = FALSE, col.names = FALSE)
write.table(tpm, file="tpm_wo_header.tab", sep="\t", na="0", quote = FALSE, col.names = FALSE)

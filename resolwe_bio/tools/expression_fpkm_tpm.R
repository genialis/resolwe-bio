#!/usr/bin/Rscript
suppressPackageStartupMessages(require('rtracklayer'))
suppressPackageStartupMessages(require('argparse'))

parser = ArgumentParser(description='Calculate expression values in FPKM and TPM units.')
parser$add_argument('annotation', help='Gene model (GTF/GFF) annotation file')
parser$add_argument('format', help='Annotation version/format')
parser$add_argument('feature_type', help='Feature type (e.g. exon)')
parser$add_argument('id_attribute', help='GFF attribute to be used as feature ID (e.g. gene_id)')
parser$add_argument('counts_file', help='Read-count data')

args = parser$parse_args(commandArgs(trailingOnly=TRUE))

# retrieve exon lengths
annot <- import.gff(args$annotation, format=args$format, genome="NA", feature.type=args$feature_type)
exon_lengths <- width(annot)
names(exon_lengths) <- elementMetadata(annot)[, args$id_attribute]

# calculate gene or transcript lengths
gene_lengths <- sapply(split(exon_lengths, names(exon_lengths)), FUN = sum)

# normalise counts by library size
ncounts <- read.delim(args$counts_file, header = FALSE)
names(ncounts) <- c("Gene", "raw_count")
ncounts$norm_count <- ncounts$raw_count / sum(ncounts$raw_count) # library_size

# obtain CPM
cpm <- ncounts$norm_count * 1e6
names(cpm) <- ncounts$Gene
cpm <- cpm[order(names(cpm))] # lexicographically order cpm

# select common genes
common_genes <- intersect(ncounts$Gene, names(gene_lengths))
subset_ncounts <- ncounts[as.character(ncounts$Gene) %in% common_genes, ]

# prepare gene_lengths in format for merging with ncounts
gene_lengths <- data.frame(Gene = names(gene_lengths), gene_lengths = unname(gene_lengths))
subset_ncounts <- merge(subset_ncounts, gene_lengths)

# normalized count normalized by gene length
subset_ncounts$norm <- with(subset_ncounts, norm_count/gene_lengths)

# obtain RPKMs
rpkm <- subset_ncounts$norm * 1e9
names(rpkm) <- subset_ncounts$Gene

# obtain TPM
tpm = exp(log(rpkm) - log(sum(rpkm)) + log(1e6))

# write data tables
write.table(cpm, file="cpm_wo_header.tab", sep="\t", na="0", quote = FALSE, col.names = FALSE)
write.table(rpkm, file="fpkm_wo_header.tab", sep="\t", na="0", quote = FALSE, col.names = FALSE)
write.table(tpm, file="tpm_wo_header.tab", sep="\t", na="0", quote = FALSE, col.names = FALSE)

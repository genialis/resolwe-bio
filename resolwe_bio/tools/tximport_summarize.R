#!/usr/bin/env Rscript
library(argparse)
library(GenomicFeatures)
library(tximport)

parser = ArgumentParser(description='Summarize TPM values to gene-level expressions.')
parser$add_argument('infile', help="Salmon output file (transcript-level).")
parser$add_argument('gtf', help='GTF file used for transcript-to-gene mapping.')
parser$add_argument('outfile', help='Output file holding the gene-level abuncance estimates.')
parser$add_argument('gene_mapping', help='Output file with transcript-to-gene mapping information.')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

# Create transcript-to-gene mapping table
db_from_gtf = makeTxDbFromGFF(args$gtf)
k = keys(db_from_gtf, keytype = "TXNAME")
tx2gene = select(db_from_gtf, k, "GENEID", "TXNAME")
# Save transcript-to-gene mapping table to file
write.table(tx2gene, file=args$gene_mapping, sep="\t", quote = FALSE, row.names = FALSE)

# Run tximport
txi_tpm = tximport(args$infile, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Ensure that the order by gene_ids is ascending,
# e.g. ENSG001, ENSG002, ..., ENSG100
tpm = data.frame(txi_tpm$abundance)
tpm = cbind(Gene = rownames(tpm), tpm)
rownames(tpm) = NULL
names(tpm) = c("Gene", "Expression")

# Write to table
write.table(tpm, file=args$outfile, sep="\t", quote = FALSE, row.names = FALSE)

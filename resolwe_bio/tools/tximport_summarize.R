#!/usr/bin/env Rscript
library(argparse)
library(GenomicFeatures)
library(tximport)

parser <- ArgumentParser(description = "Summarize TPM values to gene-level expressions.")
parser$add_argument("infile", help = "Salmon output file (transcript-level).")
parser$add_argument("gtf", help = "GTF file used for transcript-to-gene mapping.")
parser$add_argument("outfile_abundance", help = "Output file holding the gene-level abundance estimates.")
parser$add_argument("outfile_counts", help = "Output file holding the gene-level counts estimates.")
parser$add_argument("gene_mapping", help = "Output file with transcript-to-gene mapping information.")
parser$add_argument("--variance", action = "store_true", help = "Return variance in abundance estimates.")
parser$add_argument("--ignoreTxVersion", action = "store_true", help = "Strip feature_id version.")
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Create transcript-to-gene mapping table
db_from_gtf <- makeTxDbFromGFF(args$gtf)
k <- keys(db_from_gtf, keytype = "TXNAME")
tx2gene <- select(db_from_gtf, k, "GENEID", "TXNAME")
if (args$ignoreTxVersion) {
  # Strip feature_id version, e.g. ENSG000.1 -> ENSG000
  tx2gene$TXNAME <- sapply(strsplit(tx2gene$TXNAME, split = ".", fixed = TRUE), function(x) (x[1]))
  tx2gene$GENEID <- sapply(strsplit(tx2gene$GENEID, split = ".", fixed = TRUE), function(x) (x[1]))
}
# Save transcript-to-gene mapping table to file
write.table(tx2gene, file = args$gene_mapping, sep = "\t", quote = FALSE, row.names = FALSE)

# Run tximport
if (args$variance) {
  txi_tpm <- tximport(args$infile,
    type = "salmon", tx2gene = tx2gene, varReduce = TRUE,
    ignoreTxVersion = args$ignoreTxVersion
  )

  tpm_variance <- data.frame(Variance = txi_tpm$variance)
  tpm_variance <- cbind(FEATURE_ID = rownames(tpm_variance), tpm_variance)
  write.table(tpm_variance,
    file = paste0("variance_", args$outfile_abundance),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  txi_tpm <- tximport(args$infile, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = args$ignoreTxVersion)
}

# Ensure that the order by gene_ids is ascending,
# e.g. ENSG001, ENSG002, ..., ENSG100
tpm <- data.frame(Expression = txi_tpm$abundance)
tpm <- cbind(Gene = rownames(tpm), tpm)
rownames(tpm) <- NULL

counts <- data.frame(Expression = txi_tpm$counts)
counts <- cbind(Gene = rownames(counts), counts)
rownames(counts) <- NULL

# Write to table
write.table(tpm, file = args$outfile_abundance, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(counts, file = args$outfile_counts, sep = "\t", quote = FALSE, row.names = FALSE)

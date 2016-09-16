#!/usr/bin/Rscript
require('argparse') : require("DESeq2")

parser = ArgumentParser(description='Run DESeq2 differential expression.')
parser$add_argument('counts', help="Counts matrix")
parser$add_argument('--sampleConditions', nargs='+', help='design')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

conditions = unlist(strsplit(args$sampleConditions, "[,]"))
countData = read.delim(args$counts, sep="\t", header=TRUE, check.names=FALSE, row.names=1)

expDesign = data.frame(row.names = colnames(countData), condition = conditions)

dds = DESeqDataSetFromMatrix(countData = countData, colData = expDesign, design = ~ condition)
colData(dds)$condition = factor(colData(dds)$condition, levels=c("control","case"))

dds = DESeq(dds)
result = results(dds)
result = result[order(result$padj),]

head(result)

write.table(result, file='diffexp_deseq2.tab', sep='\t', quote=F, col.names=NA)

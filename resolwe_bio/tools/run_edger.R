#!/usr/bin/Rscript
library('argparse')
library("edgeR")

parser = ArgumentParser(description='Run edgeR differential expression.')
parser$add_argument('counts', help="Counts matrix")
parser$add_argument('--sampleConditions', nargs='+', help='design')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

conditions = unlist(strsplit(args$sampleConditions, "[,]"))
countData = read.delim(args$counts, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
expDesign = data.frame(row.names = colnames(countData), conditions)
group <- factor(expDesign$condition, levels=c("control","case"))

design <- model.matrix(~group, data=expDesign)
y <- DGEList(counts=countData, group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
result = topTags(lrt, n = nrow(countData))
result <- result$table
result <- result[order(result$FDR),]

write.table(result, file='diffexp_edgeR.tab', sep='\t', quote=F, col.names=NA)

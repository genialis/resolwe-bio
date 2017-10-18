#!/usr/bin/Rscript
require('argparse')
require('DESeq2')
require('tximport')

parser = ArgumentParser(description='Run DESeq2 differential expression')
parser$add_argument('--cases', nargs='+', help='Cases', required=TRUE)
parser$add_argument('--controls', nargs='+', help='Controls', required=TRUE)
parser$add_argument('--min-count-sum', type='integer', help='Minimum count sum', default=0)
parser$add_argument('--format', choices=c('rc', 'rsem'), default='rc')
args = parser$parse_args()

files <- c(args$controls, args$cases)
names(files) <- paste0('sample', 1:length(files))

conditions <- c(rep('control', length(args$controls)), rep('case', length(args$cases)))
sampleTable <- data.frame(condition=factor(conditions, levels=c('control', 'case')))

if (args$format == 'rsem') {
    txi <- tximport(files, type='rsem')
    rownames(sampleTable) <- colnames(txi$counts)
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
} else {
    data <- lapply(files, read.csv, sep='\t', row.names='Gene')
    for(i in 1:length(files)) {
        colnames(data[[i]]) <- names(files)[i]
    }
    counts <- Reduce(function(...) merge(..., by=0), data)
    cts <- counts[,-1]
    rownames(cts) <- counts[,1]
    rownames(sampleTable) <- names(files)
    dds <- DESeqDataSetFromMatrix(cts, sampleTable, ~condition)
}

dds <- dds[rowSums(counts(dds)) >= args$min_count_sum, ]
dds = DESeq(dds)

result = results(dds)
result = result[order(result$padj), ]
write.table(result, file='diffexp_deseq2.tab', sep='\t', quote=FALSE, col.names=NA)

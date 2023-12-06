#!/usr/bin/Rscript
require('argparse')
require('DESeq2')
require('tximport')

#' Raise process error
#' @param msg Error message
error <- function(msg) {
    msg <- gsub('\\n', '', msg)  # remove newline characters
    system(sprintf("re-error \"%s\"", msg))
    quit(save = "no", status = 1)
}

parser = ArgumentParser(description='Run DESeq2 differential expression')
parser$add_argument('--cases', nargs='+', help='Cases', required=TRUE)
parser$add_argument('--controls', nargs='+', help='Controls', required=TRUE)
parser$add_argument('--beta-prior', action='store_true')
parser$add_argument('--min-count-sum', type='integer', help='Minimum count sum', default=0)
parser$add_argument('--cooks-cutoff', type='double', help="Cook's cut-off", default=FALSE)
parser$add_argument('--independent', help="Independent filtering", action="store_true")
parser$add_argument('--alpha', type='double', help="Alpha", default=0.1)
parser$add_argument('--format', choices=c('rc', 'rsem', 'salmon', 'nanostring'), default='rc')
parser$add_argument('--tx2gene', help='Transcript to gene mapping')
args = parser$parse_args()

files <- c(args$controls, args$cases)
file_names <- basename(files)
names(files) <- paste0('sample', 1:length(files))

conditions <- c(rep('control', length(args$controls)), rep('case', length(args$cases)))
sampleTable <- data.frame(condition=factor(conditions, levels=c('control', 'case')))

if (args$format == 'rsem') {
    txi <- tximport(files, type='rsem')
    rownames(sampleTable) <- colnames(txi$counts)
    dds <- tryCatch(DESeqDataSetFromTximport(txi, sampleTable, ~condition), error=error)
} else if (args$format == 'salmon') {
    tx2gene = read.table(args$tx2gene)
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE)
    dds <- tryCatch(DESeqDataSetFromTximport(txi, sampleTable, ~condition), error=error)
} else if (args$format == 'nanostring'){
    data <- lapply(files, read.csv, sep='\t')
    for(i in 1:length(files)) {
        colnames(data[[i]])[2] <- names(files)[i]
    }
    counts <- Reduce(function(...) merge(..., by=1), data)
    # Round counts as the input file contains normalized counts which are doubles
    cts <- round(counts[,-1], digits = 0)
    rownames(cts) <- counts[,1]
    rownames(sampleTable) <- names(files)
    dds <- tryCatch(DESeqDataSetFromMatrix(cts, sampleTable, ~condition), error=error)
} else {
    data <- lapply(files, read.csv, sep='\t')
    for(i in 1:length(files)) {
        colnames(data[[i]])[2] <- names(files)[i]
    }
    counts <- Reduce(function(...) merge(..., by=1), data)
    cts <- counts[,-1]
    rownames(cts) <- counts[,1]
    rownames(sampleTable) <- names(files)
    dds <- tryCatch(DESeqDataSetFromMatrix(cts, sampleTable, ~condition), error=error)
}

dds <- dds[rowSums(counts(dds)) >= args$min_count_sum, ]
dds <- tryCatch(DESeq(dds, betaPrior=args$beta_prior), error=error)

result <- results(dds, cooksCutoff=args$cooks_cutoff, alpha=args$alpha,
    independentFiltering=args$independent
)
result <- result[order(result$padj), ]
write.table(result, file='diffexp_deseq2.tab', sep='\t', quote=FALSE, col.names=NA)
# Use file names for the column headers in the count matrix output file
colnames(dds) <- file_names
write.table(counts(dds), file='count_matrix.tab', sep='\t', quote=FALSE, col.names=NA)
write.table(counts(dds, normalized = TRUE), file = 'count_matrix_normalized.tab', sep ='\t', quote = FALSE, col.names = NA)

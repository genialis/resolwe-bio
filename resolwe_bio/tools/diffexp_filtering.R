#!/usr/bin/Rscript
library('argparse')

parser = ArgumentParser(description='Run edgeR differential expression.')
parser$add_argument('-counts', help="Counts matrix")
parser$add_argument('-filter', help= "Cut of value for number of reads to be filtered out")
parser$add_argument('-out', help= "Output file")
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

countData = read.delim(args$counts, sep="\t", header=TRUE, check.names=FALSE, row.names=1)

countData <- countData[apply(countData, 1, function(x) all(x >= args$filter)), ]

write.table(countData, file=args$out, sep='\t', quote=F, col.names=NA)

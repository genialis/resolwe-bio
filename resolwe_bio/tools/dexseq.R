#!/usr/bin/Rscript

require('argparse')
library("loadSubread")

parser = ArgumentParser(description='Run DEXSeq differential expression')
parser$add_argument('--cases', nargs='+', help='Cases', required=TRUE)
parser$add_argument('--controls', nargs='+', help='Controls', required=TRUE)
parser$add_argument('--expressions', help='FeatureCounts expressions', required=TRUE)
parser$add_argument('--annotation', help='DEXSeq annotation', required=TRUE)

args = parser$parse_args()

samp <- data.frame(
    row.names = c(args$controls, args$cases),
    condition = c(rep('control', length(args$controls)), rep('case', length(args$cases)))
)

dxd <- DEXSeqDataSetFromFeatureCounts(
    args$expressions,
    flattenedfile = args$annotation,
    sampleData = samp
)

dxr = DEXSeq(dxd)

write.table(dxr, file='diffexp_dexseq.tab', sep='\t', quote=FALSE, col.names=NA)
saveRDS(dxr, file='dxr.RData')

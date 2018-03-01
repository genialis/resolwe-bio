#!/usr/bin/Rscript

library('argparse')
library("loadSubread")

parser = ArgumentParser(description='Run DEXSeq differential expression')
parser$add_argument('--dxr', help='DEXSeq resuts R object', required=TRUE)
parser$add_argument('--gene', help='Gene name', required=TRUE)
parser$add_argument('-t', help='Plot transcript models', action="store_true")
parser$add_argument('-s', help='Plot individual samples', action="store_true")
parser$add_argument('-e', help='Plot exon expression', action="store_true")

args = parser$parse_args()

dxr = readRDS(args$dxr)

if (args$t) {
    if (args$s) {
        plotDEXSeq(dxr, args$gene, legend=TRUE, expression=FALSE, norCounts=TRUE, displayTranscripts=TRUE)
    }
    if (args$e) {
        plotDEXSeq(dxr, args$gene, legend=TRUE, expression=FALSE, splicing=TRUE, displayTranscripts=TRUE)
    }
    else {
        plotDEXSeq(dxr, args$gene, legend=TRUE, displayTranscripts=TRUE)
    }
} else {
    if (args$s) {
        plotDEXSeq(dxr, args$gene, legend=TRUE, expression=FALSE, norCounts=TRUE)
    }
    if (args$e) {
        plotDEXSeq(dxr, args$gene, legend=TRUE, expression=FALSE, splicing=TRUE)
    }
    else {
        plotDEXSeq(dxr, args$gene, legend=TRUE)
    }
}

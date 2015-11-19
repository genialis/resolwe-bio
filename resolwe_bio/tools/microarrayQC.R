#!/usr/bin/Rscript

require('argparse')
require('affy')
require('oligo')
require('arrayQualityMetrics')

parser = ArgumentParser(description='Run Array Quality Metrics.')
parser$add_argument('library', help='R library to use for data import.')
parser$add_argument("--logtransform", action="store_true", default=FALSE,
    help="Do logtransform of the data before running QC")
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

if (args$library == 'affy') {
    exp_data = ReadAffy()
} else {
    celFiles = list.celfiles('.', full.names = TRUE)
    exp_data = read.celfiles(celFiles)
}


if (args$logtransform) {
    arrayQualityMetrics(expressionset = exp_data,
                    outdir = "report",
                    do.logtransform = TRUE,
                    reporttitle = "QC report")
} else {
        arrayQualityMetrics(expressionset = exp_data,
                    outdir = "report",
                    do.logtransform = FALSE,
                    reporttitle = "QC report")
}

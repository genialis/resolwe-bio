#!/usr/bin/Rscript
require("Biobase")
require("limma")
require('argparse')

parser = ArgumentParser(description='Run limma differential expression.')
parser$add_argument('exp', help="Counts matrix")
parser$add_argument('--case', help='Case samples')
parser$add_argument('--control', help='Control samples')
parser$add_argument('--out', help='Output file name')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

expr = as.matrix(read.table(args$exp, header=TRUE, sep="\t", row.names=1, na.strings="#NUM!", comment.char=""))
eset <- ExpressionSet(assayData=expr)

control = unlist(strsplit(args$control, split=","))
case = unlist(strsplit(args$case, split=","))
design = apply(cbind(CONTROL=control,CASE=case), 2, as.numeric)
control
case
design

fit = lmFit(eset, design)

cont.matrix <- makeContrasts(CASEvsCONTROL=CASE-CONTROL, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output = topTable(fit2, adjust="BH", number=Inf, sort.by="P")
write.table(output, file=args$out, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

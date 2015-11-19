#!/usr/bin/Rscript
require('argparse')

parser = ArgumentParser(description='Parse closest features results file (ChIP-Seq)')
parser$add_argument('closest_features', help='Bedops closest features results file')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

parse_closest_features <- function(results_file){
	data = read.table(results_file, sep='\t', header=F, fill=T)
	split_upgeneinfo = strsplit(as.character(data$V15), '\\|')
	split_dogeneinfo = strsplit(as.character(data$V20), '\\|')
	data$upgene_orientation = sapply(split_upgeneinfo,'[', 1)
	data$dogene_orientation = sapply(split_dogeneinfo,'[', 1)
	data$upgene_distance = as.numeric(sapply(split_upgeneinfo,'[', 2))
	data$dogene_distance = as.numeric(sapply(split_dogeneinfo,'[', 2))
	data$summit_location = (1 + as.numeric(sapply(strsplit(as.character(data$V10),'\\|'),'[', 1)) + as.numeric(data$V2))
	data_filt = data[((data$upgene_orientation == '-' | data$dogene_orientation == '+') & (data$upgene_distance > 0 | data$dogene_distance > 0 )),]
	data_filt$pvalue = 10^(-data_filt$V8)
	data_filt$FDR = 10^(-data_filt$V9)
	data_filt$fold_enrichment = data_filt$V7
	data_filt$upgene_distance = data_filt$summit_location - data_filt$V12
	data_filt$upgene_distance[data_filt$upgene_orientation != '-'] = NA
	data_filt$dogene_distance = data_filt$V16 - data_filt$summit_location
	data_filt$dogene_distance[data_filt$dogene_orientation != '+'] = NA
	data_filt = data_filt[(data_filt$upgene_distance > 0 | is.na(data_filt$upgene_distance)),]
	data_filt$intergenic_distance = abs(data_filt$V12 - data_filt$V16)
	data_filt$upgene_peakScore = (data_filt$fold_enrichment)*exp(-abs(data_filt$upgene_distance/data_filt$intergenic_distance))
	data_filt$dogene_peakScore = (data_filt$fold_enrichment)*exp(-abs(data_filt$dogene_distance/data_filt$intergenic_distance))
	colnames(data_filt)[1:3] = c('chr','peakStart','peakEnd')
	data_filt$upgene_ddb_g = data_filt$V13
	data_filt$dogene_ddb_g = data_filt$V18
	out_file = paste(results_file,'_genomicContext', sep='')
	write.table(data_filt[, grep('chr|up|do|summit|peak|intergenic|pvalue|FDR|fold_enrichment', colnames(data_filt))], out_file, sep='\t', quote=F, row.names=F)
	val_ret = data_filt[, grep('chr|up|do|summit|peak|intergenic|pvalue|FDR|fold_enrichment', colnames(data_filt))]
}

parse_closest_features(args$closest_features)

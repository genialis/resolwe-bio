#!/usr/bin/Rscript
#### user inputs
#1. data
#2. replicate structure (which experiments in group A and which in group B)
### comparison always binary => two models DE: A different_from B ; NDE: A similar_to B

### mapability file: should be pre-computed and available and should be chosen based on read length of associated experiments
require('argparse') ; require('baySeq') ; require('snow')

parser = ArgumentParser(description='Run baySeq differential expression.')

parser$add_argument('expression', help='Expression file')
parser$add_argument('mappability', help='Mappability file')
parser$add_argument('--replicates', type='integer', nargs='+', help='Replicate indices')
parser$add_argument('--groups', type='integer', nargs='+', help='Group indices')

args = parser$parse_args(commandArgs(trailingOnly=TRUE))

sprintf('Expression file: %s', args$expression)
sprintf('Mappability file: %s', args$mappability)
sprintf('Replicates: %s', toString(args$replicates))
sprintf('Groups: %s', toString(args$groups))
print


getbaySeq <- function(data_set, gene_ids, replicates_, groups_, effective_gene_len) {
  print(data_set[1:3,])
  print(replicates_)
  print(groups_)
  lib_s = colSums(data_set)
  tmp_cd = new("countData", data=data_set, replicates=replicates_, libsizes=lib_s,
               groups=groups_, seglens=effective_gene_len)
  tmp_cd@annotation = gene_ids
  cl = makeCluster(4, "SOCK")
  tmp_pri = getPriors.NB(tmp_cd, samplesize=dim(data_set)[1], estimation="QL", cl=cl)
  tmp_pos = getLikelihoods.NB(tmp_pri, prs=c(0.5, 0.5), cl=cl, bootStraps=3, conv=1e-5)
  stopCluster(cl)
  val_ret = tmp_pos
}


print('Reading mappability file...')
dd_map50 = read.table(args$mappability, sep='\t', header=T)
dd_map50$effective_length = round(dd_map50$gene_len*dd_map50$mapability)

# remove genes with 0 effective length
dd_map50_filt = dd_map50[dd_map50$effective_length > 0,] ; colnames(dd_map50_filt)[1] = 'Gene'

#### read input file (raw data; matrix; all counts integers)
print('Reading gene expression data...')
inp = read.table(args$expression, sep='\t', header=T)

#### remove genes with 0 mapability
print('Removing genes with 0 mapability...')
inp_filt = merge(dd_map50_filt, inp)

## subset of data for mock run ONLY
# print('Generating subset for mock run...')
# inp_filt = inp_filt[sample((dim(inp_filt)[1]), 1000),]

#### parameters for baySeq
print('Computing baySeq parameters...')
effective_gene_len = inp_filt$effective_length
gene_ids = data.frame(inp_filt$Gene)
data_set = as.matrix(inp_filt[,grep(paste(colnames(inp)[-grep('Gene|gene|ddb_g', colnames(inp))], collapse='|'), colnames(inp_filt))])

#### replicates_ and groups_ are directly from user input
## in this example, two are 4 hour experiments and two are 0 hour experiments in this order [4, 0, 4, 0]
replicates_ = args$replicates
## groups are 'NDE: no differential expression' and 'DE: differential expression'
groups_ = list(NDE=rep(1, length(args$groups)), DE=args$groups)

#### call function to DO baySeq analysis

print('Runing baySeq...')
bays_out = getbaySeq(data_set, gene_ids, replicates_, groups_, effective_gene_len)

#### post-processing and final output

#### data relevant to the 'NDE' model
## in topCounts function, the parameter 'number' controls the number of genes that are returned
tmp_NDE = topCounts(bays_out, group="NDE", number=dim(data_set)[1]) ; colnames(tmp_NDE)[1] = 'Gene'
colnames(tmp_NDE)[grep('FDR', colnames(tmp_NDE))] = 'FDR.NDE'
colnames(tmp_NDE)[grep('Like', colnames(tmp_NDE))] = 'Lik.NDE'

#### data relevant to the 'DE' mode
tmp_DE = topCounts(bays_out, group="DE", number=dim(data_set)[1]) ; colnames(tmp_DE)[1] = 'Gene'
colnames(tmp_DE)[grep('FDR', colnames(tmp_DE))] = 'FDR.DE'
colnames(tmp_DE)[grep('Like', colnames(tmp_DE))] = 'Lik.DE'

#### final output matrix
bays_final = merge(tmp_DE, tmp_NDE)

#### write output to tab-delimited file
write.table(bays_final, 'diffexp_bayseq_partial.tab', sep='\t', quote=F, row.names=F)

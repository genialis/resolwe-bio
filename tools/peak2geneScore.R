#!/usr/bin/Rscript
require('argparse')

parser = ArgumentParser(description='Parse closest features results file (ChIP-Seq)')
parser$add_argument('peakScorefile', help='PeakScore file')
parser$add_argument('--fdr', type='double', help='FDR threshold')
parser$add_argument('--pval', type='double', help='Pval threshold')
parser$add_argument('--logratio', type='double', help='Log-ratio threshold')
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

peak2geneScore <- function(peakScorefile, fdr_th=args$fdr, pval_th=args$pval, logratio_th=args$logratio){
  d = read.table(peakScorefile, sep='\t', header=T)
  d_filt = d[unique(intersect(grep('DDB_G', d$upgene_ddb_g), grep('DDB_G', d$dogene_ddb_g))),]
  gene_id = unique(c(as.character(d_filt$upgene_ddb_g), as.character(d_filt$dogene_ddb_g)))
  gene_scores = c()
  for(gid in gene_id){
    tmp_u = d_filt[grep(gid, d_filt$upgene_ddb_g),]
    tmp_u = tmp_u[
      (tmp_u$FDR <= fdr_th & tmp_u$pvalue <= pval_th & tmp_u$fold_enrichment >= logratio_th),]
    if (dim(tmp_u)[1] > 0) { g_su = sum(tmp_u$upgene_peakScore) }
    else { g_su = 0 }
    tmp_u = d_filt[grep(gid, d_filt$dogene_ddb_g),]
    tmp_u = tmp_u[
      (tmp_u$FDR <= fdr_th & tmp_u$pvalue <= pval_th & tmp_u$fold_enrichment >= logratio_th),]
    if (dim(tmp_u)[1] > 0) { g_sd = sum(tmp_u$dogene_peakScore) }
    else { g_sd = 0 }
    g_su[is.na(g_su)] = 0 ; g_sd[is.na(g_sd)] = 0
    gene_scores = c(gene_scores, (g_su + g_sd))
  }
  val_ret = data.frame(ddb_g = gene_id, bindingScore = gene_scores)
  write.table(val_ret, "geneScore.xls", sep='\t', quote=F, row.names=F)
}

peak2geneScore(args$peakScorefile)

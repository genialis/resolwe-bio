#!/usr/bin/Rscript

library(argparse)
library(sesame)
library(GenomicRanges)

ROOT <- "/methyl_data"

#' Fetch cached manifest data
#'
#' This is needed because the process should be able to run in an
#' isolated environment without internet access. The files need to be
#' cached into the folder specified in the rootdir variable in the
#' function.
#'
#' @param x Character. Possible values are EPIC, HM27, HM450
#' @returns A list of length two where each element is a data.frame.
#'
#' $ordering: "Probe_ID" "M" "U" "DESIGN" "COLOR_CHANNEL" "col" "mask"
#' $controls: "Address" "Type" "Color_Channel" "Name"
getManifestData <- function(x) {
  rootdir <- file.path(ROOT, "platform_address")

  platform <- gsub("^(.*)\\.((hg19)|(hg38))\\.manifest$", replacement = "\\1", x = x)
  rdafile <- paste0(platform, ".address.rda")
  rdafile <- file.path(rootdir, rdafile)

  env <- new.env()
  load(rdafile, envir = env)

  adr <- paste0(platform, ".address")
  out <- get(adr, envir = env)

  return(out)
}

setExperimentHubOption("CACHE", file.path(ROOT, "ExperimentHub"))

parser <- ArgumentParser(
  description = "Parse methylation array arguments for SeSAMe pipeline."
)

parser$add_argument("-p", "--platform", help = "Platform used to process the data.", type = "character")
parser$add_argument("-m", "--manifest", help = "Which manifest file to use in importing the data", type = "character")
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

platform <- args$platform
manifest <- args$manifest

message(sprintf(fmt = "Importing manifest and control data.frame for %s from an RDS file.", manifest))
maco <- getManifestData(x = manifest)  # maco = MAnifest & COntrol

message("Importing IDAT data.")
xy <- list.files(path = "./data", full.names = TRUE)
root <- gsub("^(.*)_((Red)|(Grn)).*$", replacement = "\\1", x = xy)
root <- unique(root)

idat <- readIDATpair(
  prefix.path = root,
  platform = platform,
  manifest = maco$ordering,
  control = maco$control,
  verbose = TRUE
)

message("Running openSesame for sigset and betas.")
sigset <- openSesame(x = idat, platform = platform, what = "sigset", BPPARAM = MulticoreParam(workers = 8))
betas <- openSesame(x = idat, platform = platform, what = "beta", BPPARAM = MulticoreParam(workers = 8))

# Coerce QC data into a data frame and write to file.
message("Performing QC on the data.")
qc.attr <- c("sample_name", "mean_beta_cg", "frac_meth_cg", "frac_unmeth_cg", "GCT", "ethnicity", "sex", "age")
qc.info <- sesameQC(sset = sigset)
qc.info <- as.data.frame(qc.info)
qc.info$sample_name <- basename(root)
qc.data <- qc.info[, qc.attr]

# Rounding needed because in some cases, this is not identical between
# computers.
numcols <- sapply(qc.data, FUN = is.numeric)
qc.data[numcols] <- round(qc.data[numcols], digits = 10)

write.table(qc.data, file = "QC_data.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

message("Generating result of probe ids, betas, m-values and p-values.")
ann.betas <- as.data.frame(betas)
ann.betas$probe_ids <- rownames(ann.betas)
rownames(ann.betas) <- NULL

message("Add various variables (chr, start, end..) to probe_ids and betas.")
data.manifest <- sesameDataGet(manifest)
hgnc <- mcols(data.manifest)[, "gene_HGNC", drop = FALSE]
hgnc$probe_ids <- rownames(hgnc)
hgnc$chr <- as.character(seqnames(data.manifest))
hgnc$start <- start(data.manifest)
hgnc$end <- end(data.manifest)
hgnc$strand <- strand(data.manifest)
ann.betas <- merge(x = ann.betas, y = hgnc, all.x = TRUE)

message("Preparing beta values.")
mvals <- BetaValueToMValue(betas)
mvals <- as.data.frame(mvals)
mvals$probe_ids <- rownames(mvals)
rownames(mvals) <- NULL
# This object has two columns
#       mvals  probe_ids
# 1 -3.292382 cg00000029
# 2  4.393300 cg00000108
# 3  2.219531 cg00000109

message("Preparing p-values.")
pvals <- pval(sigset)
pvals <- as.data.frame(pvals)
pvals$probe_ids <- rownames(pvals)
rownames(pvals) <- NULL
# This object has two columns
#          pvals  probe_ids
# 1 9.611646e-04 cg01686861
# 2 1.372737e-01 cg09499248
# 3 4.319841e-05 cg00762468

message("Performing left join on betas and p-values.")
ann <- merge(x = ann.betas, y = mvals, all.x = TRUE)
ann <- merge(x = ann, y = pvals, all.x = TRUE)
colsel <- c("probe_ids", "gene_HGNC", "chr", "start", "end", "strand", "betas", "mvals", "pvals")
ann <- ann[, colsel]

message("Writing results to a gz file.")
write.table(
  x = ann,
  file = gzfile("beta_values_annotated.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
#  probe_ids         chr     start       end   strand   gene_HGNC     mvals       pvals
# <character> <character> <integer> <integer> <factor> <character> <numeric>   <numeric>
# cg00000029       chr16  53434200  53434201        -        RBL2 -3.292382    0.008586
# cg00000108        chr3  37417715  37417716        -     C3orf35  4.393300    0.000577
# cg00000109        chr3 172198247 172198248        -      FNDC3B  2.219531    0.009237
# cg00000165        chr1  90729117  90729118        +          NA -0.476389    0.032399

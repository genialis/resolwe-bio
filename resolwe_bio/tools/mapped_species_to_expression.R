#!/usr/bin/Rscript
# This script converts a mappings file (of e.g. shRNA) to a fully-fledged expression file to be used in processes
# which accept `data:expression` data. This script is part of the functionality provided by shRNAde::summaryTable()
# which works on a batch of files.

# Input parameters to this script are:
#  * file (see table (1) below on what the input looks like)
#  * sample name
#  * count_readlengths
#  * count_alignscores

# (1) count, shRNA-gene, sequence, bowtie2 alignment score
#       1 AKAP1_3352-AKAP1	GTATATTAAGTGGTGAGAATCTGGTGCC	AS:i:-12
#       1 ALKBH3_1830-ALKBH3	GTATGTAGTCTCCATTCTCTTCTGGTGCC	AS:i:0
#       1 ASF1A_5-ASF1A	GTATTATCTTCCCAATTAATGTGGATGCC	AS:i:0
#       1 ASF1B_3051-ASF1B	GTATATAAGTGATCTTTGACTGCATGCC	AS:i:-10
#       1 ASH1L_616-ASH1L	GTATTCTTTTTCCGTTTTCGTTTCTTGCC	AS:i:0

library(argparse)

default.name <- "didnt_enter"
parser <- ArgumentParser(description = "Convert a mapped species file to proper expression table.")

parser$add_argument("--input", help = "Relative or absolute path to file to be converted to expression.",
                    type = "character")
parser$add_argument("--name", help = "Name of the sample/object must be specified.",
                    type = "character",
                    default = default.name)
parser$add_argument("--readlengths",
                    help = "Species with length below this threshold will be excluded from the final count matrix.",
                    type = "integer")
parser$add_argument("--alignscores",
                    help = "Values below this threshold will be excluded from the final count matrix.",
                    type = "integer")
parser$add_argument("--output",
                    help = "Name of the output file. If not specified, defaults to`<sample_name>_count_matrix.txt`.",
                    type = "character")

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

# Do some basic sanity checks.
stopifnot(file.exists(args$input))
stopifnot(args$name != "didnt_enter")

if (is.null(args$output)) {
    output_name <- sprintf("%s_count_matrix.txt", args$name)
} else {
    output_name <- args$output
}

# Import data and put it in proper form.
xy <- read.table(args$input, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(xy) <- c("align_counts", "mapped_shRNAs", "sequences", "align_scores")
xy$read_lengths <- nchar(xy$sequences)
xy$align_scores <- as.numeric(gsub("^AS:i:", "", xy$align_scores))
xy$sample <- args$name  # append sample name

# Keep only those shRNAs which satisfy the criteria specified with arguments alignscores and readlengths.
# If no value is specified, perform no filtering.
if (is.null(args$alignscores)) args$alignscores <- min(xy$align_scores, na.rm = TRUE)
if (is.null(args$readlengths)) args$readlengths <- min(xy$read_lengths, na.rm = TRUE)

xy <- xy[xy$align_scores >= args$alignscores & xy$read_lengths >= args$readlengths, ]

# Aggreggate align_counts (sum) according to sample and shRNA.
xy <- aggregate(align_counts ~ mapped_shRNAs + sample, data = xy, FUN = sum)
  #      mapped_shRNAs    sample align_counts
  # 1 AANAT_2131-AANAT TB4445-17           30
  # 2 AANAT_2740-AANAT TB4445-17            5
  # 3  AANAT_304-AANAT TB4445-17           47
  # 4 AANAT_3349-AANAT TB4445-17          244
  # 5 AANAT_3958-AANAT TB4445-17          114
  # 6 AANAT_4567-AANAT TB4445-17           12

# Reshape data to wide format so that rows are shRNAs and columns are samples.
xy <- reshape(data = xy, timevar = "sample", idvar = "mapped_shRNAs", direction = "wide")
names(xy) <- gsub("align_counts\\.", "", names(xy))
  # mapped_shRNAs TB4445-17
  # 1 AANAT_2131-AANAT        30
  # 2 AANAT_2740-AANAT         5
  # 3  AANAT_304-AANAT        47
  # 4 AANAT_3349-AANAT       244
  # 5 AANAT_3958-AANAT       114
  # 6 AANAT_4567-AANAT        12

# Write data to file and return a list of calculated data.
write.table(xy, file = output_name, sep = "\t", row.names = FALSE, quote = FALSE)
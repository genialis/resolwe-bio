#! /bin/bash

source /etc/profile.d/resolwe-base.sh

bam_file="$1"
number_of_cores="$2"
bin_size="$3"

NAME=`basename "${bam_file}" .bam`

bamCoverage \
  --binSize "${bin_size}" \
  --bam "${bam_file}" \
  --outFileName "${NAME}.bw" \
  --outFileFormat "bigwig" \
  --numberOfProcessors "${number_of_cores}"
re-checkrc "Bam to BigWig failed."

re-save-file bigwig "${NAME}.bw"

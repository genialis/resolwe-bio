#! /bin/bash

source /etc/profile.d/resolwe-base.sh

bam_file="$1"
number_of_cores="$2"
bin_size="$3"
time_limit="$4"
species="$5"

BAM_SEQUENCES=`samtools view -c "${bam_file}"`
if [ "${BAM_SEQUENCES}" == 0 ]; then
  re-error "Bam file has no entries. No bigWig file will be made."
  exit 0
fi

NAME=`basename "${bam_file}" .bam`

time_deeptools=`python -c "import os; print(max(os.path.getsize('${bam_file}') / 1024**3, 1) * ${time_limit})"`

timeout "${time_deeptools}" \
  bamCoverage \
    --binSize "${bin_size}" \
    --bam "${bam_file}" \
    --outFileName "tmp_${NAME}.bw" \
    --outFileFormat "bigwig" \
    --numberOfProcessors "${number_of_cores}"

exit_status="$?"

if [ ${exit_status} == 124 ]
then
  re-warning "BigWig took too long to process."
elif [ ${exit_status} == 0 ]
then
  bigwig_chroms_to_ucsc.py \
    --infile "tmp_${NAME}.bw" \
    --outfile "${NAME}.bw" \
    --species "${species}"
  re-checkrc "Mapping of BigWig failed."
  re-save-file bigwig "${NAME}.bw"
else
  re-checkrc "Calculation of BigWig failed."
fi

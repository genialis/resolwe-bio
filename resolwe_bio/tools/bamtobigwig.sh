#! /bin/bash

source /etc/profile.d/resolwe-base.sh

bam_file="$1"
number_of_cores="$2"
bin_size="$3"
time_limit="$4"

NAME=`basename "${bam_file}" .bam`

time_deeptools=`python -c "import os; print(max(os.path.getsize('${bam_file}') / 1024**3, 1) * ${time_limit})"`

timeout "${time_deeptools}" \
  bamCoverage \
    --binSize "${bin_size}" \
    --bam "${bam_file}" \
    --outFileName "${NAME}.bw" \
    --outFileFormat "bigwig" \
    --numberOfProcessors "${number_of_cores}"

exit_status="$?"

if [ ${exit_status} == 124 ]
then
  re-warning "BigWig took too long to process."
elif [ ${exit_status} == 0 ]
then
  re-save-file bigwig "${NAME}.bw"
else
  re-checkrc "Calculation of BigWig failed."
fi

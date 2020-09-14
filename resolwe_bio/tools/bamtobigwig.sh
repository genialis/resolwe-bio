#! /bin/bash

source /etc/profile.d/resolwe-base.sh

bam_file="$1"
species="$2"
number_of_cores="${3:-1}"
tool="${4:-"deeptools"}"
bin_size="${5:-50}"
time_limit="${6:-480}"
normalize="${7:-None}"  # bamCoverage argument, path to .bed file
smooth="${8:-None}"  # bamCoverage argument, integer

NAME=`basename "${bam_file}" .bam`

params=(
--binSize "${bin_size}"
--bam "${bam_file}"
--outFileName "tmp_${NAME}.bw"
--outFileFormat "bigwig"
--numberOfProcessors "${number_of_cores}"
)

if [[ ${normalize} != "None" && ${tool} == "deeptools" ]]
then params+=(--normalizeUsing "${normalize}")
fi

if [[ ${smooth} != "None" && ${tool} == "deeptools" ]]
then params+=(--smoothLength "${smooth}")
fi

BAM_SEQUENCES=`samtools view -c "${bam_file}"`
if [ "${BAM_SEQUENCES}" == 0 ]; then
  re-warning "Bam file has no entries. No bigWig file will be made."
  exit 0
fi

# if BAM was created from globin reference sequence or rRNA sequence BiqWig file is not needed
samtools idxstats "${bam_file}" | cut -f 1 | head -n -1 > bam_chrom_names.txt
BAM_SOURCE=$(check_bam_source.py --infile bam_chrom_names.txt --species "${species}")
re-checkrc "Checking bam source failed."

if [ "${BAM_SOURCE}" == "globin_or_rrna" ]
then
  exit 0
fi

time=`python -c "import os; print(max(os.path.getsize('${bam_file}') / 1024**3, 1) * ${time_limit})"`

if [ "${tool}" == "deeptools" ]
then
  timeout "$((${time}/${number_of_cores}))" \
     bamCoverage "${params[@]}"

  exit_status="$?"

  if [ "${exit_status}" == 124 ]
  then
    re-warning "BigWig took too long to process."
    exit 0
  elif [ "${exit_status}" == 1 ]
  then
    re-warning "Calculation of BigWig failed."
    exit 0
  fi

else
  samtools idxstats "${bam_file}" | cut -f -2 | head -n -1 > chrom.sizes
  re-checkrc "Creating a file with chromosome sizes failed."
  timeout "$((${time}/${number_of_cores}))" \
    genomeCoverageBed \
      -bg \
      -ibam "${bam_file}" \
      -g chrom.sizes > tmp.bedgraph

    exit_status="$?"

  if [ ${exit_status} == 124 ]
  then
    re-warning "BigWig took too long to process."
    exit 0
  elif [ ${exit_status} == 1 ]
  then
    re-error "Calculation of BigWig failed."
  else
    LC_COLLATE=C sort -k1,1 -k2,2n tmp.bedgraph > tmp_sorted.bedgraph
    bedGraphToBigWig tmp_sorted.bedgraph chrom.sizes "tmp_${NAME}.bw"
  fi
fi

bigwig_chroms_to_ucsc.py \
  --infile "tmp_${NAME}.bw" \
  --outfile "${NAME}.bw" \
  --species "${species}"

re-checkrc "Mapping of BigWig failed."

re-save-file bigwig "${NAME}.bw"

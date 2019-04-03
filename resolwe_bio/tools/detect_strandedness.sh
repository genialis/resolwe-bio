#! /bin/bash

source /etc/profile.d/resolwe-base.sh

declare -A strand_codes
strand_codes=(
    ['IU']=0
    ['U']=0
    ['ISF']=1
    ['OSF']=1
    ['SF']=1
    ['ISR']=2
    ['OSR']=2
    ['SR']=2
)

BAM=$1
N_READS=$2
SALMON_INDEX=$3
THREADS=$4
FIELD=$5

if [[ $# -lt 4 || $# -gt 5 ]]; then
    re-error "Usage: detect_strandedness.sh aligned_reads n_reads salmon_index threads [field]"
fi

TOTAL_READS=$(samtools view -c "${BAM}")
SAMPLING_RATE=$(python3 -c "f = ${N_READS}/${TOTAL_READS}; print(f) if f < 1.0 else print(1.0)")

samtools view \
    -@ "${THREADS}" \
    -h \
    -s "${SAMPLING_RATE}" \
    "${BAM}" \
    | samtools sort \
    -@ "${THREADS}" \
    -n - > subsampled_sorted.bam
re-checkrc "Subsampling of alignment file (.bam) failed"

# Check if aligned reads are single or paired-end
if [[ $(samtools view -c -f 1 "${BAM}") == 0 ]]
then
    bedtools bamtofastq -i subsampled_sorted.bam -fq mate1.fastq 2> log.txt
    re-checkrc "Extraction of single-end reads from the subsampled alignment file failed"
    READS_INPUT="-r mate1.fastq"
else
    bedtools bamtofastq -i subsampled_sorted.bam -fq mate1.fastq -fq2 mate2.fastq 2> log.txt
    re-checkrc "Extraction of paired-end reads from the subsampled alignment file failed"
    READS_INPUT="-1 mate1.fastq -2 mate2.fastq"
fi

salmon --no-version-check quant \
    -i "${SALMON_INDEX}" \
    -l A \
    ${READS_INPUT} \
    -o results \
    -p "${THREADS}" \
    --minAssignedFrags 1
re-checkrc "Automated detection of strandedness failed"

if [[ -n "${FIELD}" ]]; then
    re-save-file "${FIELD}" results/lib_format_counts.json
fi

RUN="import json; print(json.load(open('results/lib_format_counts.json')).get('expected_format', 'not_detected'))"
STRAND_CODE=`python -c "${RUN}"`

if test "${strand_codes[${STRAND_CODE}]}"
    then
        STRANDEDNESS=${strand_codes[${STRAND_CODE}]}
    else
        re-error "Automated detection of strandedness failed. Re-run analysis with user-selected strandedness mode."
fi

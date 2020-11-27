#! /bin/bash

source /etc/profile.d/resolwe-base.sh

if [[ $# -lt 1 || $# -gt 5 ]]; then
  re-error "Usage: fastqc.sh reads.fastq.gz [directory [fastqc_archive [fastqc_url [fastqc-options]]]]."
fi

DIRECTORY="${2:-fastqc}"
if [ ! -d "$DIRECTORY" ]; then
  mkdir "$DIRECTORY"
fi
fastqc "$1" --extract --outdir="$DIRECTORY" "${@:5}" >&2
re-checkrc "Processing with FastQC failed."

NAME=`basename "$1" .fastq.gz`

if [ "$3" != "" ]; then
  re-save-file-list "$3" "$DIRECTORY/${NAME}_fastqc.zip"
fi

if [ "$4" != "" ]; then
  re-save-file-list "$4" "${DIRECTORY}/${NAME}_fastqc/fastqc_report.html:fastqc/${NAME}_fastqc"
fi

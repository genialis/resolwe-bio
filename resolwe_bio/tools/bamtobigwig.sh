#! /bin/bash

source /etc/profile.d/resolwe-base.sh

NAME=`basename "$1" .bam`

samtools idxstats "$1" | cut -f -2 | head -n -1 > chrom.sizes
re-checkrc
genomeCoverageBed -bg -ibam "$1" > tmp.bedgraph
re-checkrc "Bam to Bedgraph failed."
bedGraphToBigWig tmp.bedgraph chrom.sizes "${NAME}.bw"
re-checkrc "BedGraph to BigWig failed."

re-save-file bigwig "${NAME}.bw"

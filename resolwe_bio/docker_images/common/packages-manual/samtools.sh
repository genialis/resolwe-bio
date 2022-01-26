#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    samtools \
    1.14 \
    9341dabaa98b0ea7d60fd47e42af25df43a7d3d64d8e654cdf852974546b7d74 \
    https://github.com/samtools/samtools/releases/download/\${version}/samtools-\${version}.tar.bz2 \
    samtools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(samtools)

add_binary_path \
    samtools \
    samtools

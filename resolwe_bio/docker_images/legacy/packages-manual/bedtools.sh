#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    arq5x \
    bedtools \
    2.26.0 \
    65f32f32cbf1b91ba42854b40c604aa6a16c7d3b3ec110d6acf438eb22df0a4a \
    https://github.com/arq5x/bedtools2/releases/download/v\${version}/bedtools-\${version}.tar.gz \
    bedtools2

make

shopt -s extglob
rm -rf !(bin)

add_binary_path \
    arq5x \
    bedtools \
    bin

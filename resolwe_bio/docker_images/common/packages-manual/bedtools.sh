#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    quinlan-lab \
    bedtools2 \
    2.30.0 \
    333ad1ffcdc6e36005b4d6c9290677986ee97871cff92ed821c1b643d38150b8 \
    https://github.com/arq5x/bedtools2/releases/download/v\${version}/bedtools-\${version}.tar.gz \
    bedtools2

make

rm -r data
rm -r genomes
rm -r docs
rm -r tutorial

add_binary_path \
    quinlan-lab \
    bedtools2 \
    bin

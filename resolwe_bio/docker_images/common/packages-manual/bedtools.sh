#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    quinlan-lab \
    bedtools2 \
    2.29.2 \
    e3f1bf9e58740e60c3913390fe95b0c7f8fd99ceade8a406e28620448a997054 \
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

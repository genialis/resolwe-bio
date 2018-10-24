#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    quinlan-lab \
    bedtools2 \
    2.27.1 \
    c8c2acbaf3f9cb92dcf8e5cd59af6b31ae9c4598efb786ba6c84f66ca72fafd9 \
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

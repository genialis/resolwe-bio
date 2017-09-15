#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    samtools \
    1.5 \
    8542da26832ee08c1978713f5f6188ff750635b50d8ab126a0c7bb2ac1ae2df6 \
    https://github.com/samtools/samtools/releases/download/\${version}/samtools-\${version}.tar.bz2 \
    samtools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(samtools)

add_binary_path \
    samtools \
    samtools

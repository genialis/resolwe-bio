#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    samtools \
    1.10 \
    7b9ec5f05d61ec17bd9a82927e45d8ef37f813f79eb03fe06c88377f1bd03585 \
    https://github.com/samtools/samtools/releases/download/\${version}/samtools-\${version}.tar.bz2 \
    samtools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(samtools)

add_binary_path \
    samtools \
    samtools

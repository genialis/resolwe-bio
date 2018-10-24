#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    samtools \
    1.9 \
    083f688d7070082411c72c27372104ed472ed7a620591d06f928e653ebc23482 \
    https://github.com/samtools/samtools/releases/download/\${version}/samtools-\${version}.tar.bz2 \
    samtools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(samtools)

add_binary_path \
    samtools \
    samtools

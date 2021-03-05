#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    samtools \
    1.11 \
    e283cebd6c1c49f0cf8a3ca4fa56e1d651496b4d2e42f80ab75991a9ece4e5b6 \
    https://github.com/samtools/samtools/releases/download/\${version}/samtools-\${version}.tar.bz2 \
    samtools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(samtools)

add_binary_path \
    samtools \
    samtools

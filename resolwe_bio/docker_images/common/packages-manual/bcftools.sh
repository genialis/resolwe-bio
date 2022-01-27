#!/bin/bash -e

# Install bcftools from the GitHub repository

. /var/cache/build/packages-manual/common.sh

download_and_verify \
  samtools \
  bcftools \
  1.14 \
  b7ef88ae89fcb55658c5bea2e8cb8e756b055e13860036d6be13756782aa19cb \
  https://github.com/samtools/bcftools/releases/download/\${version}/bcftools-\${version}.tar.bz2 \
  bcftools-\${version}

./configure --without-curses
make

shopt -s extglob
rm -rf !(bcftools)

add_binary_path \
  samtools \
  bcftools
#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    samtools \
    bcftools \
    1.3.1 \
    12c37a4054cbf1980223e2b3a80a7fdb3fd850324a4ba6832e38fdba91f1b924 \
    https://github.com/samtools/bcftools/releases/download/\${version}/bcftools-\${version}.tar.bz2 \
    bcftools-\${version}

make pluginpath=/opt/samtools/bcftools/plugins

shopt -s extglob
rm -rf !(bcftools|vcfutils.pl|plugins)
cd plugins
rm -rf !(*.so)

add_binary_path \
    samtools \
    bcftools

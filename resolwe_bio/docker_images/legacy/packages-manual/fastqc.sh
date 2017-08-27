#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    babraham \
    fastqc \
    0.11.5 \
    dd7a5ad80ceed2588cf6d6ffe35e0f161c0d9977ed08355f5e4d9473282cbd66 \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v\${version}.zip \
    FastQC

chmod +x fastqc

add_binary_path \
    babraham \
    fastqc

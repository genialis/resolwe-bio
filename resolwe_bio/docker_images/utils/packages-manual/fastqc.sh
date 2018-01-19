#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    babraham \
    fastqc \
    0.11.7 \
    59cf50876bbe5f363442eb989e43ae3eaab8d932c49e8cff2c1a1898dd721112 \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v\${version}.zip \
    FastQC

chmod +x fastqc

add_binary_path \
    babraham \
    fastqc

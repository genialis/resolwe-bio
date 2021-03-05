#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    babraham \
    fastqc \
    0.11.9 \
    15510a176ef798e40325b717cac556509fb218268cfdb9a35ea6776498321369 \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v\${version}.zip \
    FastQC

chmod +x fastqc

add_binary_path \
    babraham \
    fastqc

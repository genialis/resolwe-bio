#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    hisat2 \
    2.2.1 \
    b7a6a4bdffabaab0ef625d047598c353b6943c1edf7c6c3a367274e0416f0bb5 \
    https://github.com/DaehwanKimLab/hisat2/archive/v\${version}.zip \
    hisat2-\${version}

# build from source
make

add_binary_path \
    ccb-jhu \
    hisat2

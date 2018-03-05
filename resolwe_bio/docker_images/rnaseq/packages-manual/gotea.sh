#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    gotea \
    0.0.4 \
    ce00cf04e1f29da153604c02ff7028a93962bf2bdbaea00755125e02012dbee8 \
    https://github.com/genialis/gotea/archive/\${version}.tar.gz \
    gotea-\${version}

make

add_binary_path \
    genialis \
    gotea

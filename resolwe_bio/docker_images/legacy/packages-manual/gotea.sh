#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    gotea \
    0.0.3 \
    eab8370b630225c215650db671c650ea662184e18825eb180c3912792bd34ab7 \
    https://github.com/genialis/gotea/archive/\${version}.tar.gz \
    gotea-\${version}

make

add_binary_path \
    genialis \
    gotea

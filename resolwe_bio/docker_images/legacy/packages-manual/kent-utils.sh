#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    encode-dcc \
    kentutils \
    302.1.0 \
    d7019c8b66c84e4c66432676c1beb799451bb389df1d3b1be5b6a8adaf22283c \
    https://github.com/ENCODE-DCC/kentUtils/archive/v\${version}.tar.gz \
    kentUtils-\${version}

rm -r src
find bin/linux.x86_64 -type f -not -name 'bedGraphToBigWig' -print0 | xargs -0 rm --

add_binary_path \
    encode-dcc \
    kentutils \
    bin/linux.x86_64

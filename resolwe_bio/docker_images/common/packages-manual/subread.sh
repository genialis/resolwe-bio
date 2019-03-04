#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    wehi \
    subread \
    1.6.3 \
    0b23094a291582725d455dbcbd3fc337e60cec04199202815b990f526b9bc8f0 \
    https://downloads.sourceforge.net/project/subread/subread-\${version}/subread-\${version}-Linux-x86_64.tar.gz \
    subread-\${version}-Linux-x86_64

add_binary_path \
    wehi \
    subread \
    bin

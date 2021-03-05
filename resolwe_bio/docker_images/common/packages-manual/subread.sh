#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    wehi \
    subread \
    2.0.1 \
    bd4d3834ea4690d27c8b06095662343ebdabd3c24c7e6c09611bfeb61ef041fa \
    https://downloads.sourceforge.net/project/subread/subread-\${version}/subread-\${version}-Linux-x86_64.tar.gz \
    subread-\${version}-Linux-x86_64

add_binary_path \
    wehi \
    subread \
    bin

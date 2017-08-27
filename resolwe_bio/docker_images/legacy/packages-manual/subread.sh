#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    wehi \
    subread \
    1.5.1 \
    008ac7e973b0d27663674b3929b05498ddb55ccd50cbb27507538588af4e1df7 \
    https://downloads.sourceforge.net/project/subread/subread-\${version}/subread-\${version}-Linux-x86_64.tar.gz \
    subread-\${version}-Linux-x86_64

add_binary_path \
    wehi \
    subread \
    bin

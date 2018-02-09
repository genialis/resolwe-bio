#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    wehi \
    subread \
    1.6.0 \
    62b505910803ad4bdc54d77da14b201c909fd90d27e7a94d89d1297b13f7541b \
    https://downloads.sourceforge.net/project/subread/subread-\${version}/subread-\${version}-Linux-x86_64.tar.gz \
    subread-\${version}-Linux-x86_64

add_binary_path \
    wehi \
    subread \
    bin

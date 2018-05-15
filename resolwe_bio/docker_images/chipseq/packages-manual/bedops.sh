#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bedops \
    bedops \
    2.4.32 \
    cce96ded8e0276ed02a1b08f949c0ced09991305eedb61f5e2aa6f79e6e80b1f \
    https://github.com/bedops/bedops/releases/download/v\${version}/bedops_linux_x86_64-v\${version}.tar.bz2

add_binary_path \
    bedops \
    bedops \
    bin

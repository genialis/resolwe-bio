#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bedops \
    bedops \
    2.4.15 \
    d18171695bfa41525185f1ea4ed82a20920b2dce8984fea0b6b5e7dbde722fda \
    https://github.com/bedops/bedops/releases/download/v\${version}/bedops_linux_x86_64-v\${version}.tar.bz2

add_binary_path \
    bedops \
    bedops \
    bin

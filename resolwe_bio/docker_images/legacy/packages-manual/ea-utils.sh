#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ea \
    ea-utils \
    1.04.807 \
    aa09d25e6aa7ae71d2ce4198a98e58d563f151f8ff248e4602fa437f12b8d05f \
    https://github.com/ExpressionAnalysis/ea-utils/archive/\${version}.tar.gz \
    ea-utils-\${version}/clipper

make

add_binary_path \
    ea \
    ea-utils

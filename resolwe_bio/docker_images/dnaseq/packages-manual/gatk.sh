#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    gatk \
    4.0.11.0 \
    5ee23159be7c65051335ac155444c6a49c4d8e3515d4227646c0686819934536 \
    https://github.com/broadinstitute/gatk/releases/download/\${version}/gatk-\${version}.zip \
    gatk-\${version}

add_binary_path \
    broadinstitute \
    gatk

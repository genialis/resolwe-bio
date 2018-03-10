#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    gatk \
    4.0.2.1 \
    b8f6f2d39a5400b00f6430b4474df0ce3c36b4ddfd10e6e7a03537e87e9c3b76 \
    https://github.com/broadinstitute/gatk/releases/download/\${version}/gatk-\${version}.zip \
    gatk-\${version}

add_binary_path \
    broadinstitute \
    gatk

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    gatk \
    4.1.2.0 \
    ffc5f9b3d4b35772ee5dac3060b59dc657f30e830745160671d84d732c30dc65 \
    https://github.com/broadinstitute/gatk/releases/download/\${version}/gatk-\${version}.zip \
    gatk-\${version}

add_binary_path \
    broadinstitute \
    gatk

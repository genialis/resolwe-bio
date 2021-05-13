#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    gatk \
    4.2.0.0 \
    dd11cc8e3bc7a23c2c226366428f0908c902765eabbc1c641c736c06b80aaf78 \
    https://github.com/broadinstitute/gatk/releases/download/\${version}/gatk-\${version}.zip \
    gatk-\${version}

add_binary_path \
    broadinstitute \
    gatk

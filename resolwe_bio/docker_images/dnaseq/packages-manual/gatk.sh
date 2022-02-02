#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    gatk \
    4.2.4.1 \
    2dc007f534cca1e29253864601e480ce24ce376b4b2bbd7a491571f9ecedc2e2 \
    https://github.com/broadinstitute/gatk/releases/download/\${version}/gatk-\${version}.zip \
    gatk-\${version}

add_binary_path \
    broadinstitute \
    gatk

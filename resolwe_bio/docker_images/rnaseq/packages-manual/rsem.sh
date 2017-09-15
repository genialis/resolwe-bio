#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    deweylab \
    rsem \
    1.3.0 \
    a6776aae25a690375504974d65e8554d6e3ddab93cef5b82969064854d22c09f \
    https://github.com/deweylab/RSEM/archive/v\${version}.tar.gz \
    RSEM-\${version}

make

add_binary_path \
    deweylab \
    rsem

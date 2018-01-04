#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    deweylab \
    rsem \
    1.3.0 \
    ecfbb79c23973e1c4134f05201f4bd89b0caf0ce4ae1ffd7c4ddc329ed4e05d2 \
    https://github.com/deweylab/RSEM/archive/v\${version}.tar.gz \
    RSEM-\${version}

make

add_binary_path \
    deweylab \
    rsem

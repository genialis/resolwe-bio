#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    sunnyisgalaxy \
    moabs \
    1.3.2 \
    a9a29d7625757954b38d41cb14b048169d6026bcbf932ff9b16148af30040887 \
    https://github.com/sunnyisgalaxy/moabs/archive/\${version}.zip \
    moabs-\${version}

find . -depth -not \( -path './bin*' -o \( -type d -not -empty \) \) -delete

add_binary_path \
    sunnyisgalaxy \
    moabs \
    bin

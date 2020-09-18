#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kreuger \
    trimgalore \
    0.6.6 \
    b8db8ffd131d9d9e7c8532a5a1f1caee656c0c58d3eafd460fee3c39b9fcab5e \
    https://github.com/FelixKrueger/TrimGalore/archive/\${version}.tar.gz \
    TrimGalore-\${version}

add_binary_path \
    kreuger \
    trimgalore

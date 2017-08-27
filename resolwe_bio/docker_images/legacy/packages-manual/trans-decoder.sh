#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    transdecoder \
    transdecoder \
    3.0.1 \
    753a5fac5bfd04466aeabff48053c92e876cece8906b96de3b72f23f86fafae7 \
    https://github.com/TransDecoder/TransDecoder/archive/v\${version}.tar.gz \
    TransDecoder-\${version}

make
rm -r sample_data

add_binary_path \
    transdecoder \
    transdecoder

add_binary_path \
    transdecoder \
    transdecoder \
    util

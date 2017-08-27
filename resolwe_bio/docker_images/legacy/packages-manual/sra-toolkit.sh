#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    nih \
    sratoolkit \
    2.8.0 \
    ce229dd9eaebdc7d257a1e70566c95eb970aac2017c86beb51162df8b275b54b \
    ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/\${version}/sratoolkit.\${version}-ubuntu64.tar.gz \
    sratoolkit.\${version}-ubuntu64

add_binary_path \
    nih \
    sratoolkit \
    bin

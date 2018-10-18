#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    smithlab \
    methpipe \
    3.4.3 \
    56716370211a7b45b0a3a2994afb64d64c15dd362028f9ecd8a0551a6e6d65c3 \
    http://smithlabresearch.org/downloads/methpipe-\${version}.tar.bz2 \
    methpipe-\${version}

make
make install

shopt -s extglob
rm -rf !(bin)

add_binary_path \
    smithlab \
    methpipe \
    bin

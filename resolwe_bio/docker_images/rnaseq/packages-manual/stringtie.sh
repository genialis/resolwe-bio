#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    stringtie \
    2.0 \
    00229938ee771ec5cd3596e846a83fc9f82aebd08cc5e6287f5fa23860830413 \
    https://github.com/gpertea/stringtie/releases/download/v\${version}/stringtie-\${version}.Linux_x86_64.tar.gz \
    stringtie-\${version}.Linux_x86_64 \
    stringtie-\${version}.tar.gz

add_binary_path \
    ccb-jhu \
    stringtie

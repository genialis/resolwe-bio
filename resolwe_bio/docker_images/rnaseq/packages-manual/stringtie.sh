#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    stringtie \
    1.3.6 \
    5905874f484dc7086c18e1f7bed8164d6f38119ab45d28db05dd9134c9a22f02 \
    http://ccb.jhu.edu/software/stringtie/dl/stringtie-\${version}.Linux_x86_64.tar.gz \
    stringtie-\${version}.Linux_x86_64 \
    stringtie-\${version}.tar.gz

add_binary_path \
    ccb-jhu \
    stringtie

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    hisat2 \
    2.1.0 \
    7f68df59db7db93b8a4710d934d75577ac9138ecccc2a8b26d700925354c09bc \
    ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-\${version}-Linux_x86_64.zip \
    hisat2-\${version}

# remove debugging files, documentation and examples
rm *-debug
rm -r doc
rm -r example

add_binary_path \
    ccb-jhu \
    hisat2

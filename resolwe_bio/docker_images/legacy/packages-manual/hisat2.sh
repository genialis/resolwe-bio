#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    hisat2 \
    2.0.3-beta \
    4b32770619efcc6c14dcc9c39544cd4f7946cd286776b29965347f11f00e699d \
    ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-\${version}-Linux_x86_64.zip \
    hisat2-\${version}

# remove debugging files, documentation and examples
rm *-debug
rm -r doc
rm -r example

add_binary_path \
    ccb-jhu \
    hisat2

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kent \
    wigtobigwig \
    no_version \
    69c1af385a8a767444f59838e5f2dde12319e4005f512b126801d0860038b216 \
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig#.elf

mv wigtobigwig wigToBigWig

add_binary_path \
    kent \
    wigtobigwig

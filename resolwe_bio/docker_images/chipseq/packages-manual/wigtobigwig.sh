#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kent \
    wigtobigwig \
    no_version \
    9b3afbbdca5a0227cbc9d7d82de5f83017b04ab53524040c159fe6db8a726f71 \
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig#.elf

mv wigtobigwig wigToBigWig

add_binary_path \
    kent \
    wigtobigwig

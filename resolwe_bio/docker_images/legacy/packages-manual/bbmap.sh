#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bbmap \
    bbmap \
    36.99 \
    331e534ef4652e948ee6f8eb4cae286a35301f81f81b4060358acd9e72fde19b \
    https://downloads.sourceforge.net/project/bbmap/BBMap_\${version}.tar.gz \
    bbmap

add_binary_path \
    bbmap \
    bbmap

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bbmap \
    bbmap \
    37.90 \
    f2df007b581aec1b4278e2c115f52afef2bb3e4117c1cd5bd11bbf4993d10bbe \
    https://downloads.sourceforge.net/project/bbmap/BBMap_\${version}.tar.gz \
    bbmap

add_binary_path \
    bbmap \
    bbmap

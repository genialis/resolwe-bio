#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    prinseq \
    prinseq-lite \
    0.20.4 \
    9b5e0dce3b7f02f09e1cc7e8a2dd77c0b133e5e35529d570ee901f53ebfeb56f \
    https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-\${version}.tar.gz \
    prinseq-lite-\${version}

find -iname "*.pl" -type f | xargs chmod 0755

add_binary_path \
    prinseq \
    prinseq-lite

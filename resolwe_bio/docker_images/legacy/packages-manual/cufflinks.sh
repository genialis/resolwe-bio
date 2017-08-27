#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ctl \
    cufflinks \
    2.2.1 \
    9615bba771e6c89546644c653a99127e49c4453396f2cc7f98b6363c157bb961 \
    'https://data2.genialis.com/f/fa3cf97573b74569a5fe/?dl=1#.tar.bz2' \
    cufflinks

add_binary_path \
    ctl \
    cufflinks \
    bin

add_library_path \
    ctl \
    cufflinks \
    lib

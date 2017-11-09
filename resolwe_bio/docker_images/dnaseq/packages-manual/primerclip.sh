#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    primerclip \
    v171018 \
    ff118b659e9b40f8847645126b87c158fe04bc7c7f6f709c0ce50c1f0bb46f70 \
    'https://data2.genialis.com/f/8e8fb294e097476aad75/?dl=1#.elf'

add_binary_path \
    genialis \
    primerclip

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    primerclip \
    2 \
    5b74d96bab9bff0b7880ca134c4d67c67230e7d36f161269282d5994254da3c4 \
    'https://data2.genialis.com/f/6fde14e212bd471fa7fc/?dl=1#.elf'

add_binary_path \
    genialis \
    primerclip

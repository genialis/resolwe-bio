#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    primerclip \
    v180309 \
    9496caea64395e5159a0f03abb093a3a79f6d3f0e590ed5abf8af2dabb5682f0 \
    'https://data2.genialis.com/f/5286d798d4d64acfb5d2/?dl=1#.elf'

add_binary_path \
    genialis \
    primerclip

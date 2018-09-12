#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genialis \
    primerclip \
    v171018 \
    ff118b659e9b40f8847645126b87c158fe04bc7c7f6f709c0ce50c1f0bb46f70 \
    'https://drive.google.com/uc?id=1xlgMb3aaVUOIMVXqSZIS7ulrOJ3p3BiA#.elf'

add_binary_path \
    genialis \
    primerclip

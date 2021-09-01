#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bwa-mem2 \
    bwa-mem2 \
    2.2.1 \
    b4cfdbce8cc07cdf3f6a920facabc29c976cf77dd53573369508111d6d1c555b \
    https://github.com/bwa-mem2/bwa-mem2/releases/download/v\${version}/bwa-mem2-\${version}_x64-linux.tar.bz2 \
    bwa-mem2-\${version}_x64-linux


add_binary_path \
    bwa-mem2 \
    bwa-mem2

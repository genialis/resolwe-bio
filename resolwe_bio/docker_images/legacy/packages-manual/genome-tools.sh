#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genometools \
    genometools \
    1.5.3 \
    aa25bb55c2aeba14954dedfff2b359d22e2219ce62ee23056d10483c713d7650 \
    https://github.com/genometools/genometools/archive/v\${version}.tar.gz \
    genometools-\${version}

make 64bit=yes cairo=no -j $(nproc)
make 64bit=yes cairo=no install
cd ..
rm -rf genometools

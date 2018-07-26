#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    griffithlab \
    regtools \
    0.5.0 \
    c9a5f05b23ffb219c35d12a0403b34db7f7cee2b0be7dc3b6f71154dc838917e \
    https://github.com/griffithlab/regtools/archive/\${version}.tar.gz \
    regtools-\${version}

mkdir build
cd build/
cmake ..
make

shopt -s extglob
rm -rf !(regtools)
cd ..
rm -rf !(build)

add_binary_path \
    griffithlab \
    regtools \
    build


#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    alexdobin \
    star \
    2.5.2a \
    2a372d9bcab1dac8d35cbbed3f0ab58291e4fbe99d6c1842b094ba7449d55476 \
    https://github.com/alexdobin/STAR/archive/\${version}.tar.gz \
    STAR-\${version}

rm -r doc
rm -r source
rm -r extras
rm -r bin/Linux_x86_64_static
rm -r bin/MacOSX_x86_64

add_binary_path \
    alexdobin \
    star \
    bin/Linux_x86_64

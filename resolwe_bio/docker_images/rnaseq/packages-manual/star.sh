#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    alexdobin \
    star \
    2.5.3a \
    2a258e77cda103aa293e528f8597f25dc760cba188d0a7bc7c9452f4698e7c04 \
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

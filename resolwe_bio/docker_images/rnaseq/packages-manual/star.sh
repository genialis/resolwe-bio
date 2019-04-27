#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    alexdobin \
    star \
    2.7.0f \
    daead0a87b259e6311f48f31ccaab5bd53ef1da9d6e5cd7210e0d714b7e23204 \
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

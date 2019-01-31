#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    alexdobin \
    star \
    2.7.0a \
    d0771c620176b8d53df5ced818b4dd0f2c8ac537bc44828cdfefddc0fd05caa1 \
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

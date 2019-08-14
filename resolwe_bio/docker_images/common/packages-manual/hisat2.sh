#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ccb-jhu \
    hisat2 \
    0f01dc6397a \
    53ca877e4c9af975d58301643d5745f680ce391850b33f3c8db5a62a59f518d8 \
    https://github.com/DaehwanKimLab/hisat2/archive/\${version}.zip \
    hisat2-\${version}

# build from source
make

# remove documentation and examples
rm -r doc
rm -r example

add_binary_path \
    ccb-jhu \
    hisat2

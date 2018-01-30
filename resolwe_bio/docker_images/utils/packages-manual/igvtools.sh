#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    igvtools \
    2.3.98 \
    07027c179f25960bab9919f255c0f8e08f0861861ac6dc02d92be8313e0316a4 \
    http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_\${version}.zip \
    IGVTools

rm *.txt
rm -r genomes

add_binary_path \
    broadinstitute \
    igvtools

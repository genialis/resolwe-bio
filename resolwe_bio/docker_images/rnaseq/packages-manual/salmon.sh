#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    combine-lab \
    salmon \
    0.14.0 \
    d1c4a969e735f21050cb3a9595e581a9ff21bb94a3c744cb61a6ab8079b66750 \
    https://github.com/COMBINE-lab/salmon/releases/download/v\${version}/Salmon-\${version}_linux_x86_64.tar.gz \
    salmon-latest_linux_x86_64

rm sample_data.tgz

add_binary_path \
    combine-lab \
    salmon \
    bin

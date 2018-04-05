#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    combine-lab \
    salmon \
    0.9.1 \
    67d19cec3ad4df4f6b24f1d64be39f445227e197c451e2eeb2aa0b71ee59d093 \
    https://github.com/COMBINE-lab/salmon/releases/download/v\${version}/Salmon-\${version}_linux_x86_64.tar.gz \
    Salmon-latest_linux_x86_64

rm sample_data.tgz

add_binary_path \
    combine-lab \
    salmon \
    bin

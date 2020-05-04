#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    combine-lab \
    salmon \
    1.2.1 \
    898417a36828fc1b64456276584a12beaf19342dbd0a921134dbd0e22ac7dcb5 \
    https://github.com/COMBINE-lab/salmon/releases/download/v\${version}/Salmon-\${version}_linux_x86_64.tar.gz \
    salmon-latest_linux_x86_64

rm sample_data.tgz

add_binary_path \
    combine-lab \
    salmon \
    bin

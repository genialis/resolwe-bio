#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    nimblegen \
    hsqutils \
    1.0 \
    bb54448bb4a0edf02f6c9b071e3e299ce410dc5d6e4ca036074dac4ed4d494c7 \
    https://github.com/NimbleGen/bioinformatics/releases/download/v1.0/hsqutils_v1_0.zip \
    hsqutils_v1_0

mkdir bin
echo '#!/bin/sh' > bin/hsqutils
echo 'exec java -Xmx4g -Xms4g -jar /opt/nimblegen/hsqutils/hsqutils.jar "$@"' > bin/hsqutils
chmod +x bin/hsqutils

add_binary_path \
    nimblegen \
    hsqutils \
    bin

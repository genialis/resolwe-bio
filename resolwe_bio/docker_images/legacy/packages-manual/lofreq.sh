#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    lofreq \
    lofreq \
    2.1.2 \
    e226e72a1c67148367cde26322053cbd9381ec446129828c10bd0f9269e26623 \
    https://sourceforge.net/projects/lofreq/files/lofreq_star-\${version}_linux-x86-64.tgz \
    lofreq_star-\${version}

add_binary_path \
    lofreq \
    lofreq \
    bin

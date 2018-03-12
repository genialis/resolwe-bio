#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    lofreq \
    lofreq \
    2.1.3.1 \
    508ad72a36c31984feb204a6e4861e39157a98229cacb10cf20287c33efdd2b9 \
    https://github.com/CSB5/lofreq/raw/f2c00dc3345f9be2a5f400d171cfe942e335cf5f/dist/lofreq_star-\${version}_linux-x86-64.tgz \
    lofreq_star-\${version}

add_binary_path \
    lofreq \
    lofreq \
    bin

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    tommyau \
    bamclipper \
    1.1.1 \
    081a6008b16209ed27bb7fd6e1e86e89efc431e8cb958990f84a969b301485e2 \
    'https://github.com/tommyau/bamclipper/archive/master.zip'

mv bamclipper-master bamclipper

add_binary_path \
    tommyau \
    bamclipper \
    bamclipper

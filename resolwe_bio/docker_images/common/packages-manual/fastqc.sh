#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    babraham \
    fastqc \
    0.11.8 \
    ca87fe77807e4ac796b6cad949858921fd20652c4038f586f05ece94b5022129 \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v\${version}.zip \
    FastQC

chmod +x fastqc

add_binary_path \
    babraham \
    fastqc

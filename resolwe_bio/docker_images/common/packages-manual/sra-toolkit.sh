#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ncbi \
    sra-toolkit \
    2.10.0 \
    0a283d8cf7ff9da9ff836056bd0a7e4c49c583cdc3c440e2c316c4dbb3029c80 \
    https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/\${version}/sratoolkit.\${version}-ubuntu64.tar.gz \
    sratoolkit.\${version}-ubuntu64

# create a config file
mkdir -p /root/.ncbi
printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg
printf '/config/default = "false"\n' >> /root/.ncbi/user-settings.mkfg
printf '/repository/user/main/public/cache-disabled = "true"\n' >>/root/.ncbi/user-settings.mkfg

add_binary_path \
    ncbi \
    sra-toolkit \
    bin

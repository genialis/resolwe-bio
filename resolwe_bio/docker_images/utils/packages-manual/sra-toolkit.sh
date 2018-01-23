#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ncbi \
    sra-toolkit \
    2.8.2-1 \
    42b2fc111f55a26be2b245b79da3b57f53f0b7ab56daa71d28eb6652371476f3 \
    https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    sratoolkit.\${version}-ubuntu64

mkdir /root/.ncbi/
echo '## auto-generated configuration file - DO NOT EDIT ##' >>/root/.ncbi/user-settings.mkfg
echo '' >>/root/.ncbi/user-settings.mkfg
echo '/config/default = "false"' >>/root/.ncbi/user-settings.mkfg
echo '/repository/user/cache-disabled = "true"' >>/root/.ncbi/user-settings.mkfg
echo '/repository/user/default-path = "/root/ncbi"' >>/root/.ncbi/user-settings.mkfg

add_binary_path \
    ncbi \
    sra-toolkit \
    bin

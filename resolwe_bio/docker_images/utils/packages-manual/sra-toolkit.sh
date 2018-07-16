#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ncbi \
    sra-toolkit \
    2.9.0 \
    911301a6f85bdbf8e786a6b6561d32f262cd356ae52acd27830e39fb07e41b12 \
    https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/\${version}/sratoolkit.\${version}-ubuntu64.tar.gz \
    sratoolkit.\${version}-ubuntu64

mkdir /.ncbi/
echo '## auto-generated configuration file - DO NOT EDIT ##' >>/.ncbi/user-settings.mkfg
echo '' >>/.ncbi/user-settings.mkfg
echo '/config/default = "false"' >>/.ncbi/user-settings.mkfg
echo '/repository/user/cache-disabled = "true"' >>/.ncbi/user-settings.mkfg

add_binary_path \
    ncbi \
    sra-toolkit \
    bin

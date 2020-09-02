#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ncbi \
    sra-toolkit \
    2.10.8 \
    3718cd77147a842538ea04475a0bb3052e525fdeeb07f4440f6d2794cc477083 \
    https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/\${version}/sratoolkit.\${version}-ubuntu64.tar.gz \
    sratoolkit.\${version}-ubuntu64

# create a config file
mkdir -p /home/sra-toolkit/.ncbi
printf '/LIBS/GUID = "%s"\n' `uuidgen` > /home/sra-toolkit/.ncbi/user-settings.mkfg
printf '/config/default = "false"\n' >> /home/sra-toolkit/.ncbi/user-settings.mkfg
printf '/repository/user/main/public/cache-disabled = "true"\n' >> /home/sra-toolkit/.ncbi/user-settings.mkfg

# cp -R /root/.ncbi /opt/ncbi/

add_binary_path \
    ncbi \
    sra-toolkit \
    bin

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    ensembl \
    vep \
    104.3 \
    832c73d831f9a62349b55e21efe1b74d65498280b4d021d03f19f3acd9d7307b \
    https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/\${version}.zip \
    ensembl-vep-release-\${version}


perl INSTALL.pl --AUTO al --NO_UPDATE
# this command disables the removal of these two packages with apt-get autoremove
apt-mark manual libmysqlclient21 mysql-common

add_binary_path \
    ensembl \
    vep

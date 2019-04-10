#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

vendor=10x
name=cellranger
version=3.0.2
fileid=1HYqfowTCgZE6S3iDUR82I0gqeccORa_t

echo "Downloading '${vendor}/${name}' version '${version}'..."

mkdir -p "/opt/${vendor}/${name}"
cd "/opt/${vendor}"

# Token is needed for large files since Google Drive can't scan them
# for viruses and direct download is thus disabled
wget --quiet --save-cookies cookies.txt 'https://docs.google.com/uc?export=download&id='${fileid} -O- \
     | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt

wget --quiet --load-cookies cookies.txt -O ${fileid} \
     'https://docs.google.com/uc?export=download&id='${fileid}'&confirm='$(<confirm.txt)

tar -xf "${fileid}" --directory "${name}"

rm -f confirm.txt cookies.txt "${fileid}"

add_binary_path \
    ${vendor} \
    ${name} \
    ${name}-${version}

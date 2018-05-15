. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genometools \
    genometools \
    1.5.10 \
    a6aa7f158a3cef90fea8d0fe24bfad0c3ee96b17b3ba0c1f6462582593af679e \
    https://github.com/genometools/genometools/archive/v\${version}.tar.gz \
    genometools-\${version}

make cairo=no -j $(nproc)
make cairo=no install
cd ..
rm -rf genometools

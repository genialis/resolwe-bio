. /var/cache/build/packages-manual/common.sh

download_and_verify \
    genometools \
    genometools \
    1.6.2 \
    974825ddc42602bdce3d5fbe2b6e2726e7a35e81b532a0dc236f6e375d18adac \
    https://github.com/genometools/genometools/archive/v\${version}.tar.gz \
    genometools-\${version}

make cairo=no -j $(nproc)
make cairo=no install
cd ..
rm -rf genometools

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kundajelab \
    phantompeakqualtools \
    1.2 \
    86cbcca80b65f150b1cdbea673d8a47caba88c2db6b3b567a80f2c797c9a1668 \
    https://github.com/kundajelab/phantompeakqualtools/archive/\${version}.tar.gz \
    phantompeakqualtools-\${version}

sed -i '1i#!/usr/bin/env Rscript\nlibrary(caTools)\n' run_spp.R
mv run_spp.R spp
chmod +x spp
rm README.md spp_*

add_binary_path \
    kundajelab \
    phantompeakqualtools

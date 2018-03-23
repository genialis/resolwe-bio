#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    bokeh \
    assets \
    0.12.13 \
    205b00f79345046f6a10e6c362806e38be6e5910cf9f942a044e9198d669fa96 \
    https://cdn.pydata.org/bokeh/release/bokeh-\${version}.min.css

download_and_verify \
    bokeh \
    assets \
    0.12.13 \
    6b2ff777f56e7ecec30448a03354fbf0c7b747200f987d94ece9410160b9f99e \
    https://cdn.pydata.org/bokeh/release/bokeh-\${version}.min.js

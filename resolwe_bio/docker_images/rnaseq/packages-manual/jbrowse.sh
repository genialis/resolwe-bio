#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    gmod \
    jbrowse \
    1.16.6 \
    4e196691b80e7b84dcb1d56191e4f18bf4b7c2b543c3d5e9c58a30178837cec4 \
    https://github.com/GMOD/jbrowse/releases/download/\${version}-release/JBrowse-\${version}.zip \
    JBrowse-\${version}

# patch setup.sh script to prevent formatting of example data and building
# support for legacy tools
sed -i '/Formatting Volvox example data .../,$d' setup.sh
./setup.sh
# remove all files and directories except those we explicitly want to keep
find . -depth -not \( \
    -path './bin*' -o \
    -path './src/perl5*' -o \
    -path './extlib/lib/perl5*' \
    -o \( -type d -not -empty \) \
\) -delete
# JBrowse's Perl scripts want to create files during runtime, however, we
# don't know which UID/GID will be used, so we just give permissions to
# everyone
chmod --recursive 777 extlib

add_binary_path \
    gmod \
    jbrowse \
    bin

#!/bin/bash -e

#
# Licensing conditions stated in http://hgdownload.soe.ucsc.edu/admin/exe/
# ========================================================================
#
# (updated Mar 29, 2018)
#
# The UCSC Genome Browser and Blat software are free for academic, nonprofit, and personal use. A
# license is required for commercial download and installation of these binaries, with the exception
# of items built from the following source code directories, which are freely available for all uses:
#
#  - kent/src/utils (includes *big* tools)
#  - kent/src/lib
#  - kent/src/inc
#  - kent/src/tagStorm
#  - kent/src/parasol
#  - kent/src/htslib
#  - kent/src/hg/autoSql
#  - kent/src/hg/autoXml
#
# Note
# ====
#
# The bigBedToBed tool resides in the kent/src/utils source directory.
#


. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kent \
    bedGraphToBigWig \
    kent-v365 \
    6dc7daaf18bd7b14933eb32cbbee559f5cd69d636578caf2b29e13996f9445c1 \
    https://drive.google.com/uc?id=1GjBuK1JxSHrxWZhHnrTZmMsxwX6xGUWs#.elf

add_binary_path \
    kent \
    bedGraphToBigWig

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
    kent-v364 \
    296caadb7e60645cdc689df13f55a8bd07991dd597de75574f929d399d6cbbf7 \
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig#.elf

add_binary_path \
    kent \
    bedGraphToBigWig

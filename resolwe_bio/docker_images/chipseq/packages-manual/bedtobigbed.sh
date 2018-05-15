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
    bedToBigBed \
    kent-v364 \
    f20ec4f9cb2d740191554a68fefe86610aed998ec1cdc90b6016c3954e825846 \
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed#.elf

wget https://raw.githubusercontent.com/ENCODE-DCC/kentUtils/master/src/hg/lib/encode/narrowPeak.as

add_binary_path \
    kent \
    bedToBigBed

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
    kent-v365 \
    a66027e1a51f803fef6bf234b90d0e4d6deee6027acfe2c0a65712f3157cca5d \
    https://drive.google.com/uc?id=1N2pUaJuHgCUS4gztTcoybrRWw1m8WNtv#.elf

wget https://raw.githubusercontent.com/ENCODE-DCC/kentUtils/master/src/hg/lib/encode/narrowPeak.as

add_binary_path \
    kent \
    bedToBigBed

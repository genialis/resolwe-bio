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
# The wigToBigWig tool resides in the kent/src/utils source directory.
#

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kent \
    wigToBigWig \
    kent-v365 \
    278a7da9210b1b82c20dbce33c263094cd7f14b20fdc7f5dbe330bf709ddd04e \
    https://drive.google.com/uc?id=19M6SZmS0xxaKe3srxt01d0FrdaIYWnHe#.elf

add_binary_path \
    kent \
    wigToBigWig

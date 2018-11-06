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
# The bedClip tool resides in the kent/src/utils source directory.
#


. /var/cache/build/packages-manual/common.sh

download_and_verify \
    kent \
    bedClip \
    kent-v373 \
    fa4b081042d7013e0e6f928dde02a6b2f5d4760410811257b8508bf329fc1e55 \
    https://drive.google.com/uc?id=1_Bfhnvbs72W4UMeATPoCRFcunOIy0Ktb#.elf

add_binary_path \
    kent \
    bedClip

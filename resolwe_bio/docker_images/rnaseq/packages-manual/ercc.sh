#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# Table with information about ERCC spike-ins
# Original file is available at https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt

download_and_verify \
    resolwebio \
    assets \
    ERCC_table \
    5151b394332477301851404d37b907cb6f27c9aaa2569e288531a839d6ca2b3b \
    https://drive.google.com/uc?id=193PpxSoNhFwYUpTMh9nGt44K72tc0F9f \
    '' \
    ERCC_table.txt

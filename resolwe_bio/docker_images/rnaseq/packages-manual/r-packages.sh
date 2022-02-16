#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# Install R packages.
#d XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-r.txt | \
    Rscript --slave --no-save --no-restore-history -e "renv::install(readLines('stdin'))"

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
Rscript --slave --no-save --no-restore-history -e " \
  devtools::install_github('jkokosar/RNASeqT'); \
  devtools::install_github('jvrakor/Subread_to_DEXSeq', subdir = 'loadSubread'); \
  renv::install('genialis/shRNAde'); \
"

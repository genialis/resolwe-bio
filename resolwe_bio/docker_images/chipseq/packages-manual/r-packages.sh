#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
Rscript --slave --no-save --no-restore-history -e " \
  install.packages(c('spp', 'devtools'), dependencies=TRUE); \
  install.packages('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz', repos=NULL, type='source'); \
"

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.

Rscript --slave --no-save --no-restore-history -e " \
  devtools::install_github('genialis/ChIPQC') \
"
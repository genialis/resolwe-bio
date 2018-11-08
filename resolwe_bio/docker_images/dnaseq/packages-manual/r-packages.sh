#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
Rscript --slave --no-save --no-restore-history -e " \
  library(devtools); \
  install_github('jkokosar/chemut') \
"

#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# These steps need to be performed because certain packages produce a "database corrupt error". The solution is to
# restart R, but this is not feasible during the docker built or building of dependencies. This chunk also installs
# packages from bioconductor. This entire step should be considered temporary until we can get packrat to work for
# packages from all universes (CRAN, Bioconductor).
Rscript --slave --no-save --no-restore-history -e "install.packages('assertthat')"

Rscript --slave --no-save --no-restore-history -e \
"source('https://bioconductor.org/biocLite.R'); \
BiocInstaller::biocLite(c('IRanges', 'GenomicRanges', 'Rsamtools'))\
"

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
Rscript --slave --no-save --no-restore-history -e \
"devtools::install_github('jkokosar/chemut@436547b6cc36faad3e8c26bd382bbf383efc1b65')"

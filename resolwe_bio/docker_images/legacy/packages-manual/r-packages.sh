#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# Install R packages.
# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-r.txt | \
    Rscript --slave --no-save --no-restore-history \
        -e "install.packages(readLines('stdin'), repos='https://cloud.r-project.org')"

# XXX: This is unverifiable and in WITHOUT HTTPS and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-r-bioconductor.txt | \
    Rscript --slave --no-save --no-restore-history \
        -e "source('http://www.bioconductor.org/biocLite.R'); biocLite(readLines('stdin'))"

# XXX: This is unverifiable and thus may compromise the whole image.
# XXX: Use notary (https://github.com/ropenscilabs/notary) when ready.
Rscript --slave --no-save --no-restore-history -e " \
  library(devtools); \
  install_github('jkokosar/chemut'); \
  install_github('jkokosar/RNASeqT') \
"

download_and_verify \
    genialis \
    mirna30 \
    0 \
    cc312fd80f73c0efadd4d104d86806709db2c28382193d1e9f694a8a7abedaa7 \
    'https://data2.genialis.com/f/d12b143bb5d44210b835/?dl=1#.cdf'

Rscript --slave --no-save --no-restore-history -e " \
    library(methods); \
    library('makecdfenv'); \
    make.cdf.package('mirna30.cdf', species ='H_sapiens', packagename='mirna30cdf') \
"
R CMD INSTALL mirna30cdf
cd
rm -rf /opt/genialis/mirna30

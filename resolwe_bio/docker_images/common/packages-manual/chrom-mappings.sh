#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# All files are available on Genialis Drive and were aquired with create_mapping_file.py
# available on genialis-scripts GitHub Repository.

download_and_verify \
    chrom_mappings \
    assets \
    GRCh38.p12_ensembl2UCSC \
    1722bba20ba1565370cc7d2e0bacf3947b3da213b62efb72bf12c62b8e2df06e \
    https://drive.google.com/uc?id=1TAPSkvYfe2036WFLxBwdPMM0KM57OvBc \
    '' \
    GRCh38.p12_ensembl2UCSC.txt

download_and_verify \
    chrom_mappings \
    assets \
    GRCh38.p12_NCBI2UCSC \
    e3eaa5fa31301e8f28035ed88552b803e31a4c93a5063f8e78d864193eff77c0 \
    https://drive.google.com/uc?id=1Xn4VlyayBSJ0sx-stsF52EH4PgaM1GIQ \
    '' \
    GRCh38.p12_NCBI2UCSC.txt

download_and_verify \
    chrom_mappings \
    assets \
    GRCm38.p6_ensembl2UCSC \
    566988707c6ae198266d6bbda38caff16014186ad08d0db28809e9bc2d73ae4d \
    https://drive.google.com/uc?id=1i20f8tuqwJBAzCDYwpRjK1zfsZGaXRft \
    '' \
    GRCm38.p6_ensembl2UCSC.txt

download_and_verify \
    chrom_mappings \
    assets \
    GRCm38.p6_NCBI2UCSC \
    035eabcc5f367260122b0e34d72a0464515d7ec02dc0fb0ba329f6732e004203 \
    https://drive.google.com/uc?id=1o3LO9-C7hkygHRvAjELJXbttF5O5Z7tu \
    '' \
    GRCm38.p6_NCBI2UCSC.txt


download_and_verify \
    chrom_mappings \
    assets \
    rn6_ensembl2UCSC \
    a274bea86b73cc0d96e24e0b4506286f4fd28d4ee9590afa07836478b2ce344b \
    https://drive.google.com/uc?id=1u4AuJQfeR4pw9P6wYeFbnhu_JT2kIDkU \
    '' \
    rn6_ensembl2UCSC.txt

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

# rRNA chromosome names for Homo sapiens
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Homo_sapiens_rRNA \
    4ea6bf407616ad0e934494be51fe767755a5bc07f76e6d0ee3f30eb13c070760 \
    https://drive.google.com/uc?id=1R6vBnBh-05yoUFfYyDkCTGhfJ06yrdah \
    '' \
    Homo_sapiens_rRNA.txt

# rRNA chromosome names for Mus musculus
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Mus_musculus_rRNA \
    775760d053fdabe4248425db6bb3181fe56433151ea5719f0e376608a57635ed \
    https://drive.google.com/uc?id=1MNRvOk8nPjPXsSkZqIEFEZ9QlNLhvUes \
    '' \
    Mus_musculus_rRNA.txt

# rRNA chromosome names for Rattus_norvegicus
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Rattus_norvegicus_rRNA \
    2263632cdd7960058f5d304d29d455d502fbe311c1f69783c6f8b5cca08fa522 \
    https://drive.google.com/uc?id=1OkHsmJRnKYTk_8QTpA3W1w7qrJaOr2Ei \
    '' \
    Rattus_norvegicus_rRNA.txt

# globin chromosome names for Homo sapiens
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Homo_sapiens_globin_RNA \
    5b4a699e4ea6ee9d473f1e542ba0002929730d7e55e2db1148947d5b56fa9c7c \
    https://drive.google.com/uc?id=1dTEXjxE814egst94NB5HaO6tTom0NXxA \
    '' \
    Homo_sapiens_globin_RNA.txt

# globin chromosome names for Mus musculus
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Mus_musculus_globin_RNA \
    142bd5f98df669c15d09b47df9eb2a9c4dd7a918abda7dd2094fba9b523d68eb \
    https://drive.google.com/uc?id=1zcb0KQYPYWxVDcK_FQCvrYQsusnUyc9h \
    '' \
    Mus_musculus_globin_RNA.txt

# globin chromosome names for Rattus norvegicus
download_and_verify \
    chrom_mappings \
    globin_rrna \
    Rattus_norvegicus_globin_RNA \
    2263632cdd7960058f5d304d29d455d502fbe311c1f69783c6f8b5cca08fa522 \
    https://drive.google.com/uc?id=1NX4BmIsq72K8guYXqLnq1trfJdUzg822 \
    '' \
    Rattus_norvegicus_globin_RNA.txt

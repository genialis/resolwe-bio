FROM docker.io/resolwe/base:ubuntu-16.04

LABEL MAINTAINER Resolwe Bioinformatics authors https://github.com/genialis/resolwe-bio

RUN apt-get update && \

    apt-get -y install --no-install-recommends \
      fastqc \
      python \
      python-cutadapt \
      python-htseq \
      python-pip \
      r-base-core \
      samtools \
      unzip \
      wget \
      #Required for R packages.
      build-essential \
      libz-dev \
      libcurl4-openssl-dev \
      libxml2-dev \
      && \

    echo "Installing HISAT2..." && \
    HISAT2_VERSION=2.1.0 && \
    HISAT2_SHA256SUM=7f68df59db7db93b8a4710d934d75577ac9138ecccc2a8b26d700925354c09bc && \
    wget -q ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -O /tmp/HISAT2.zip && \
    echo "$HISAT2_SHA256SUM */tmp/HISAT2.zip" | sha256sum -c - && \
    mkdir -p /opt/johnshopkins/ && \
    unzip -q /tmp/HISAT2.zip -d /opt/johnshopkins/ && \
    rm /tmp/HISAT2.zip && \
    echo "PATH=\$PATH:/opt/johnshopkins/hisat2-$HISAT2_VERSION" >> /etc/profile.d/rnaseq.sh && \

    echo "Installing RSEM..." && \
    RSEM_VERSION=1.3.0 && \
    RSEM_SHA256SUM=a6776aae25a690375504974d65e8554d6e3ddab93cef5b82969064854d22c09f && \
    wget -q https://github.com/deweylab/RSEM/archive/v$RSEM_VERSION.tar.gz -O /tmp/RSEM.tar.gz && \
    echo "$RSEM_SHA256SUM */tmp/RSEM.tar.gz" | sha256sum -c - && \
    mkdir -p /opt/deweylab/RSEM-$RSEM_VERSION && \
    tar -xf /tmp/RSEM.tar.gz --directory /opt/deweylab/RSEM-$RSEM_VERSION --strip-components=1 && \
    rm /tmp/RSEM.tar.gz && \
    ORIG_DIR=$PWD && \
    cd /opt/deweylab/RSEM-$RSEM_VERSION && \
    make && \
    make install && \
    cd $ORIG_DIR && \
    echo "PATH=\$PATH:/opt/deweylab/RSEM-$RSEM_VERSION" >> /etc/profile.d/rnaseq.sh && \

    echo "Installing STAR..." && \
    STAR_VERSION=2.5.3a && \
    STAR_SHA256SUM=2a258e77cda103aa293e528f8597f25dc760cba188d0a7bc7c9452f4698e7c04 && \
    wget -q https://github.com/alexdobin/STAR/archive/$STAR_VERSION.tar.gz -O /tmp/STAR.tar.gz && \
    echo "$STAR_SHA256SUM */tmp/STAR.tar.gz" | sha256sum -c - && \
    mkdir -p /opt/alexdobin/STAR-$STAR_VERSION && \
    tar -xf /tmp/STAR.tar.gz --directory /opt/alexdobin/STAR-$STAR_VERSION --strip-components=1 && \
    rm /tmp/STAR.tar.gz && \
    rm -r /opt/alexdobin/STAR-$STAR_VERSION/doc && \
    rm -r /opt/alexdobin/STAR-$STAR_VERSION/source && \
    rm -r /opt/alexdobin/STAR-$STAR_VERSION/extras && \
    rm -r /opt/alexdobin/STAR-$STAR_VERSION/bin/Linux_x86_64 && \
    rm -r /opt/alexdobin/STAR-$STAR_VERSION/bin/MacOSX_x86_64 && \
    echo "PATH=\$PATH:/opt/alexdobin/STAR-$STAR_VERSION/bin/Linux_x86_64_static" >> /etc/profile.d/rnaseq.sh && \

    echo "Installing R packages..." && \
    Rscript --slave --no-save --no-restore-history -e "\
      package_list = c( \
        'argparse' \
      ); \
      install.packages(package_list, repos='https://cloud.r-project.org') \
    " && \

    echo "Installing Bioconductor R packages..." && \
    Rscript --slave --no-save --no-restore-history -e " \
      package_list = c( \
        'rtracklayer' \
      ); \
      source('http://www.bioconductor.org/biocLite.R'); \
      biocLite(package_list) \
    " && \

    # Cleanup
    apt-get autoremove --purge -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /tmp/*


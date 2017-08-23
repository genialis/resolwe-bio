FROM docker.io/resolwe/base:ubuntu-16.04

LABEL MAINTAINER Resolwe Bioinformatics authors https://github.com/genialis/resolwe-bio

RUN apt-get update && \

    apt-get -y install --no-install-recommends \
      python \
      python-pip \
      r-base-core \
      samtools \
      unzip \
      wget \
      # Required for fastqc.
      default-jre \
      perl \
      #Required for R packages.
      build-essential \
      libz-dev \
      libcurl4-openssl-dev \
      libxml2-dev \
      && \
    
    # Required for HTSeq.
    echo "Installing numpy..." && \
    pip install numpy==1.13.1 && \

    # Required for HTSeq.
    echo "Installing pysam..." && \
    pip install pysam==0.11.2.2 && \

    echo "Installing HTSeq..." && \
    pip install HTSeq==0.9.1 && \

    echo "Installing cutadapt..." && \
    pip install cutadapt==1.14 && \

    echo "Installing FastQC..." && \
    FastQC_VERSION=0.11.5 && \
    FastQC_SHA256SUM=dd7a5ad80ceed2588cf6d6ffe35e0f161c0d9977ed08355f5e4d9473282cbd66 && \
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$FastQC_VERSION.zip -O /tmp/FastQC.zip && \
    echo "$FastQC_SHA256SUM */tmp/FastQC.zip" | sha256sum -c - && \
    mkdir -p /opt/babraham/FastQC-$FastQC_VERSION && \
    unzip -q /tmp/FastQC.zip -d /opt/babraham/FastQC-$FastQC_VERSION && \
    rm /tmp/FastQC.zip && \
    chmod +x /opt/babraham/FastQC-$FastQC_VERSION/FastQC/fastqc && \
    echo "PATH=\$PATH:/opt/babraham/FastQC-$FastQC_VERSION/FastQC" >> /etc/profile.d/rnaseq.sh && \

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

    echo "Installing eXpress..." && \
    express_VERSION=1.5.1 && \
    express_SHA256SUM=50e268af69cd5f561ae238e7d155ebb4a6505219bd93bc25a7667f7bf92b7bd5 && \
    wget -q https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz -O /tmp/express.tgz && \
    echo "$express_SHA256SUM */tmp/express.tgz" | sha256sum -c - && \
    mkdir -p /opt/pachterlab/express-$express_VERSION && \
    tar -xf /tmp/express.tgz --directory /opt/pachterlab/express-$express_VERSION --strip-components=1 && \
    rm /tmp/express.tgz && \
    rm -r /opt/pachterlab/express-$express_VERSION/sample_data && \
    echo "PATH=\$PATH:/opt/pachterlab/express-$express_VERSION" >> /etc/profile.d/rnaseq.sh && \

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


FROM docker.io/resolwe/base:ubuntu-14.04

LABEL MAINTAINER Resolwe Bioinformatics authors https://github.com/genialis/resolwe-bio

RUN apt-get update && \
    apt-get -y install --no-install-recommends software-properties-common && \

    # required to install python2.6
    add-apt-repository -y ppa:fkrull/deadsnakes && \
    apt-get update && \

    apt-get -y install --no-install-recommends \
      build-essential \
      python \
      python-dev \
      python-pip \
      wget \
      samtools \
      # required for bamliquidator
      libbam-dev \
      libboost-dev \
      libboost-timer-dev \
      libgoogle-perftools-dev \
      libhdf5-7 \
      libhdf5-serial-dev \
      libtbb2 \
      libtbb-dev \
      libtcmalloc-minimal4 \
      python-numpy \
      python-scipy \
      python-pandas \
      python-redis \
      python-software-properties \
      python-tables \
      python-numexpr \
      # required for ROSE2
      r-base \
      # required for MACS1.4
      python2.6 \
      && \

    echo "Installing bokeh..." && \
    pip install bokeh==0.9.3 && \

    echo "Installing ROSE2..." && \
    pip install rose2==1.0.5 && \

    echo "Installing bamplot..." && \
    pip install bamplot==0.9.0 && \

    echo "Installing MACS2..." && \
    pip install MACS2==2.1.1.20160309 && \

    echo "Installing MACS..." && \
    MACS_VERSION=1.4.2 && \
    MACS_SHA256SUM=44032de245f22f8656be1b884342f0ad643bdd2a40568773d97f9345914f6301 && \
    wget -q https://cloud.github.com/downloads/taoliu/MACS/macs_$MACS_VERSION.deb -O /tmp/macs.deb && \
    echo "$MACS_SHA256SUM */tmp/macs.deb" | sha256sum -c - && \
    dpkg -i /tmp/macs.deb && \
    rm /tmp/macs.deb && \

    echo "Installing bamliquidator..." && \
    BAMLIQUIDATOR_VERSION=v1.3 && \
    BAMLIQUIDATOR_SHA256SUM=049094f633b99040e8d4be683bf2b1be01cafc6383987a034e5421e6822a823d && \
    wget -q https://github.com/BradnerLab/pipeline/archive/bamliquidator_$BAMLIQUIDATOR_VERSION.tar.gz -O /tmp/bamliquidator.tar.gz && \
    echo "$BAMLIQUIDATOR_SHA256SUM */tmp/bamliquidator.tar.gz" | sha256sum -c - && \
    mkdir -p /opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION && \
    tar -xf /tmp/bamliquidator.tar.gz --directory /opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION --strip-components=1 && \
    cd /opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION/bamliquidator_internal && \
    make && \
    cd / && \
    rm /tmp/bamliquidator.tar.gz && \
    echo "PATH=\$PATH:/opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION/bamliquidator_internal" \
      >> /etc/profile.d/bamliquidator.sh && \
    echo "PATH=\$PATH:/opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION/bamliquidator_internal/bamliquidatorbatch"  \
      >> /etc/profile.d/bamliquidator.sh && \
    ln -s /opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION/bamliquidator_internal/bamliquidatorbatch/bamliquidator_batch.py \
      /opt/BradnerLab/bamliquidator-$BAMLIQUIDATOR_VERSION/bamliquidator_internal/bamliquidatorbatch/bamliquidator_batch && \

    # Cleanup
    apt-get purge -y \
      build-essential \
      python-dev \
      libhdf5-serial-dev \
      libboost-dev \
      libboost-timer-dev \
      libgoogle-perftools-dev \
      libtbb-dev && \
    apt-get autoremove --purge -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /tmp/*

CMD ["/bin/bash"]

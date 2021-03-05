FROM public.ecr.aws/genialis/resolwebio/base:ubuntu-20.04-03042021

MAINTAINER Resolwe Bioinformatics authors https://github.com/genialis/resolwe-bio

# Force matplotlib backend to a headless one.
ENV MPLBACKEND Agg

# Additional files:
COPY assets/ /opt/resolwebio/assets/

# Cleanup any previus build requirement files except common.sh file
ONBUILD RUN find /var/cache/build/ ! -type d ! -name 'common.sh' -delete

# Package manifests and custom build scripts.
ONBUILD ADD packages-*.txt /var/cache/build/
ONBUILD ADD packages-manual /var/cache/build/packages-manual/

ONBUILD RUN export DEBIAN_FRONTEND=noninteractive && \
    # Install packages.
    apt-get update && \
    if [ -f /var/cache/build/packages-debian.txt ]; then \
    echo "Installing distribution packages (runtime requirements)..." && \
    sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-debian.txt | \
    xargs apt-get --no-install-recommends -y install \
    ; \
    fi && \
    if [ -f /var/cache/build/packages-debian-buildreq.txt ]; then \
    echo "Installing distribution packages (build requirements)..." && \
    sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-debian-buildreq.txt | \
    xargs apt-get --no-install-recommends -y install \
    ; \
    fi && \
    # Install Python packages.
    if [ -f /var/cache/build/packages-python2.txt ]; then \
    echo "Installing Python 2 packages..." && \
    pip install --require-hashes --requirement /var/cache/build/packages-python2.txt \
    ; \
    fi && \
    if [ -f /var/cache/build/packages-python2-stage2.txt ]; then \
    echo "Installing stage 2 Python 2 packages..." && \
    pip install --require-hashes --requirement /var/cache/build/packages-python2-stage2.txt \
    ; \
    fi && \
    if [ -f /var/cache/build/packages-python3.txt ]; then \
    echo "Installing Python 3 packages..." && \
    pip3 install --require-hashes --requirement /var/cache/build/packages-python3.txt \
    ; \
    fi && \
    if [ -f /var/cache/build/packages-python3-stage2.txt ]; then \
    echo "Installing stage 2 Python 3 packages..." && \
    pip3 install --require-hashes --requirement /var/cache/build/packages-python3-stage2.txt \
    ; \
    fi && \
    # Build manual packages.
    if [ -f /var/cache/build/packages-manual.txt ]; then \
    echo "Building manual packages..." && \
    sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-manual.txt | \
    xargs --max-args 1 -I '{}' bash -c 'bash /var/cache/build/packages-manual/{}.sh || exit 255' \
    ; \
    fi && \
    # Cleanup packages.
    echo "Cleaning up packages..." && \
    if [ -f /var/cache/build/packages-debian-buildreq.txt ]; then \
    sed -e 's/^#.*$//g' -e '/^$/d' /var/cache/build/packages-debian-buildreq.txt | \
    xargs apt-get -y purge \
    ; \
    fi && \
    apt-get autoremove --purge -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /tmp/*

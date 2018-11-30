FROM docker.io/resolwebio/common:1.0.0

MAINTAINER Resolwe Bioinformatics authors https://github.com/genialis/resolwe-bio

# Force matplotlib backend to a headless one.
ENV MPLBACKEND Agg

# Additional files:
COPY assets/ /opt/resolwebio/assets/

# Set the default JAVA version. Required for compatibility with GATK4.
RUN update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

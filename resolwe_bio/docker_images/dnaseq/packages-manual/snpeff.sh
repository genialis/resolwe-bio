#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    snpeff \
    snpeff \
    4_3t \
    d55a7389a78312947c1e7dadf5e6897b42d3c6e942e7c1b8ec68bb35d2ae2244 \
    https://sourceforge.net/projects/snpeff/files/snpEff_v\${version}_core.zip \
    snpEff

rm -rf clinEff examples
java -jar snpEff.jar download GRCh37.75 -v

mkdir bin

cat << 'EOF' >bin/snpEff
#!/bin/sh

java -Xmx16g -jar /opt/snpeff/snpeff/snpEff.jar "$@"
EOF
chmod +x bin/snpEff

cat << 'EOF' >bin/SnpSift
#!/bin/sh

java -Xmx16g -jar /opt/snpeff/snpeff/SnpSift.jar "$@"
EOF
chmod +x bin/SnpSift

# Ensure data is downloaded during build. If not, fail.
test -d data

add_binary_path \
    snpeff \
    snpeff \
    bin

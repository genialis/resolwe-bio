#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    snpeff \
    snpeff \
    5_1 \
    01a244c8043fd4abe0e8087d0d3f1f7e05bdd6edc71df7a315d55cac59fb285d \
    https://snpeff.blob.core.windows.net/versions/snpEff_v\${version}_core.zip \
    snpEff

rm -rf clinEff examples

java -jar snpEff.jar download GRCh37.75 -v
java -jar snpEff.jar download GRCh38.99 -v

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

add_binary_path \
    snpeff \
    snpeff \
    scripts

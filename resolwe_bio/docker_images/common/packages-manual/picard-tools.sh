#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    broadinstitute \
    picard-tools \
    2.26.10 \
    84a056466cf23155e2126b33fde61b5485217b365d4899c1bdc82fe9284f2395 \
    https://github.com/broadinstitute/picard/releases/download/\${version}/picard.jar

echo 'alias picard-tools="java -jar /opt/broadinstitute/picard-tools/picard.jar"' >> /etc/profile
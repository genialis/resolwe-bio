#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    hartleys \
    qorts \
    1.3.0 \
    036523a1b024121244f6a1aea5b41eb0f83d8619fb3b598fa125d2999b00fda8 \
    https://github.com/hartleys/QoRTs/releases/download/\${version}/QoRTs.jar


cat <<'EOF' >QoRTs
#!/bin/bash

java "${1}" -jar /opt/hartleys/qorts/QoRTs.jar "${@:2}"
EOF
chmod +x QoRTs

add_binary_path \
    hartleys \
    qorts

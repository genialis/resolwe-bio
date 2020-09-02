#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

download_and_verify \
    usadellab \
    trimmomatic \
    0.36 \
    4846c42347b663b9d6d3a8cef30da2aec89fc718bf291392c58e5afcea9f70fe \
    http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-\${version}.zip \
    Trimmomatic-\${version}

cat <<'EOF' >TrimmomaticSE
#!/bin/bash

java -jar /opt/usadellab/trimmomatic/trimmomatic-0.36.jar SE "$@"
EOF
chmod +x TrimmomaticSE

cat <<'EOF' >TrimmomaticPE
#!/bin/bash

java -jar /opt/usadellab/trimmomatic/trimmomatic-0.36.jar PE "$@"
EOF
chmod +x TrimmomaticPE

add_binary_path \
    usadellab \
    trimmomatic

#
# Download package and verify SHA256 digest.
#
download_and_verify() {
    local vendor=$1
    local name=$2
    local version=$3
    local sha256sum=$4
    # Eval is required so that the url and subdir arguments can contain references
    # to the \${version} variable.
    local url=$(eval "echo \"$5\"")
    local subdir=$(eval "echo \"$6\"")

    local filename=$(basename "$url")
    local extension="${filename#*.}"

    echo "Downloading '${vendor}/${name}' version '${version}'..."

    mkdir -p "/opt/${vendor}/${name}"
    cd "/opt/${vendor}"

    echo "Fetching '${url}'..."
    wget --tries=3 --retry-connrefused --timeout=30 --quiet "${url}" --output-document "${name}.${extension}"

    echo "Verifying package..."
    set +e
        echo "${sha256sum} ${name}.${extension}" | sha256sum -c --status -
        local verify_status=$?
    set -e

    if (( ${verify_status} )); then
        >&2 echo "ERROR: SHA256 digest mismatch."
        rm -f "${name}.${extension}"
        exit 1
    fi

    # Unpack.
    echo "Unpacking..."
    case "$extension" in
        *.tar*|tar*|*.tgz|tgz)
            tar -xf "${name}.${extension}" --directory ${name}
            ;;
        *.zip|zip)
            unzip -q "${name}.${extension}" -d ${name}
            ;;
        *.elf|elf)
            # The .elf extension should be used when downloading raw binaries.
            mv "${name}.${extension}" "${name}/${name}"
            chmod +x "${name}/${name}"
            ;;
        *.cdf|cdf)
            mv "${name}.${extension}" "${name}/"
            ;;
        *)
            >&2 echo "ERROR: Unknown archive format '${extension}'."
            exit 1
    esac

    rm -f "${name}.${extension}"
    cd "${name}"

    if [[ -n ${subdir} ]]; then
        mv "${subdir}"/* .
        rm -rf "${subdir}"
    fi
}

#
# Ensure binaries from the given directory are in default search path.
#
add_binary_path() {
    local vendor=$1
    local name=$2
    local bindir=$3

    echo "PATH=\$PATH:/opt/${vendor}/${name}/${bindir}" >> /etc/profile.d/resolwebio-base.sh
}

#
# Ensure libraries from the given directory are in default search path.
#
add_library_path() {
    local vendor=$1
    local name=$2
    local libdir=$3

    echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/${vendor}/${name}/${libdir}" >> /etc/profile.d/resolwebio-base.sh
}

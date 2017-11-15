#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

# TODO: Add support for specifying PPAs without the need for manual packages.
add-apt-repository ppa:bradner-computation/pipeline
apt-get update
apt-get -y install bamliquidator


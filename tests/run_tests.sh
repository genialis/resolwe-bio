#!/bin/bash
#
# A script for running Resolwe tests with Jenkins.
#
# To use this script, add an "Execute shell" "Build" step to your project and
# put in the following:
# ./tests/run_tests.sh
#

set -x

pushd $WORKSPACE

rm -rf .venv
virtualenv .venv
source .venv/bin/activate
pip install --process-dependency-links .

./tests/manage.py test resolwe_bio

cloc --exclude-dir=.venv,reports --by-file --xml --out=reports/cloc.xml .

popd

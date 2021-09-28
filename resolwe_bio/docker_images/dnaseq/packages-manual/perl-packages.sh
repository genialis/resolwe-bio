#!/bin/bash -e

. /var/cache/build/packages-manual/common.sh

perl -MCPAN -e 'install(DBI);install(DBD::mysql);install(Archive::Zip);install(Module::Build);install(Try::Tiny)'
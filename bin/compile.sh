#!/usr/bin/env bash
RETCODE=0;
./bin/build_dependencies.sh
RETCODE=$RETCODE+$?;
./bin/compile_lib.sh "$@"
RETCODE=$RETCODE+$?;
./bin/compile_bin.sh "$@"
RETCODE=$RETCODE+$?;
exit $(($RETCODE))

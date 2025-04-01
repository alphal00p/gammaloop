#!/usr/bin/env bash
RETCODE=0;
./bin/build_dependencies.sh
RETCODE=$RETCODE+$?;
if [ ! $(($RETCODE)) == 0 ]
then
    echo "ERROR: could not build dependencies";
    exit $(($RETCODE))
fi
./bin/compile_lib.sh --release "$@"
RETCODE=$RETCODE+$?;
if [ ! $(($RETCODE)) == 0 ]
then
    echo "ERROR: could not build compile gammaloop library";
    exit $(($RETCODE))
fi
./bin/compile_bin.sh --release "$@"
RETCODE=$RETCODE+$?;
if [ ! $(($RETCODE)) == 0 ]
then
    echo "ERROR: could not build compile gammaloop binary";
    exit $(($RETCODE))
fi
exit $(($RETCODE))

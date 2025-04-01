#!/usr/bin/env bash
RETCODE=0;
echo "Building pyo3 library... ["$@"]"
cargo build "$@" --lib --features="python_api";
RETCODE=$RETCODE+$?;
rm -f ./python/gammaloop/_gammaloop.so;
if [[ $1 = "--release" ]]  || [[ $1 = "--profile=release" ]]
then
    COMPILE_PATH="release"
else
    COMPILE_PATH="debug"
fi
if [ "$(uname)" == "Darwin" ]; then
    ln -s ../../target/$COMPILE_PATH/lib_gammaloop.dylib ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
else
    ln -s ../../target/$COMPILE_PATH/lib_gammaloop.so ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
fi
if [ ! -f ./python/gammaloop/_gammaloop.so ]; then
  echo -e "\n\033[0;31mERROR: Could not find target pyo3 library in rust build '$COMPILE_PATH'.\033[0m\n"
  RETCODE=$RETCODE+1
fi

exit $(($RETCODE))

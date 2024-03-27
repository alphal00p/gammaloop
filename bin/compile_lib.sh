#!/usr/bin/env bash
RETCODE=0;
echo "Building pyo3 library..."
cargo build --release "$@" --lib --features="python_api";
RETCODE=$RETCODE+$?;
rm -f ./python/gammaloop/_gammaloop.so;
if [ "$(uname)" == "Darwin" ]; then
    ln -s ../../target/release/lib_gammaloop.dylib ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
else
    ln -s ../../target/release/lib_gammaloop.so ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
fi
exit $(($RETCODE))

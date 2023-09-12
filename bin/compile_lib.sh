#!/usr/bin/env bash
cargo build $@ --lib --features="python_api";
rm -f ./python/gammaloop/_gammaloop.so;
if [ "$(uname)" == "Darwin" ]; then
    ln -s ../../target/release/lib_gammaloop.dylib ./python/gammaloop/_gammaloop.so;
else
    ln -s ../../target/release/lib_gammaloop.so ./python/gammaloop/_gammaloop.so;
fi

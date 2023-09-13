#!/usr/bin/env bash

./bin/build_dependencies.sh
./bin/compile_lib.sh $@
./bin/compile_bin.sh $@
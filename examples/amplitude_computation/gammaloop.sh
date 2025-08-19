#!/usr/bin/env bash

rm -rf gammaloop_state
rm -rf TMP
../../target/$1/cli -g $2_generation.yaml -r $2_runtime.yaml -o import model ../../src/test_resources/gammaloop_models/sm.yaml
../../target/$1/cli -g $2_generation.yaml -r $2_runtime.yaml -o import amplitude ./$3.dot
../../target/$1/cli -o generate
../../target/$1/cli inspect --process-id 0 --name $3 $4

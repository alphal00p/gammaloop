#!/usr/bin/env bash

rm -rf gammaloop_state
rm -rf TMP
../../target/$1/cli -g $2_generation.yaml -r $2_runtime.yaml -o import model ../../src/test_resources/gammaloop_models/sm.yaml
../../target/$1/cli -g $2_generation.yaml -r $2_runtime.yaml -o import amplitude ./$3.dot
../../target/$1/cli -o generate
../../target/$1/cli -o inspect --process-id 0 --name $3 $4

# ./gammaloop.sh release valentin_settings qqx_aaa_subtracted "-m -p 0.1 0.2 0.3"
# gives:
# For input point xs:
#
#( 0.005812358643936371, 0.17620819117478337, 0.9008918628686365 )
#
#The evaluation of integrand 'qqx_aaa_subtracted' is:
#
#( 1.1948744538168098e-1, -1.9693624149798648e1 i)

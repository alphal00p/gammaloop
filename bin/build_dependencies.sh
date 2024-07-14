cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b nested_evaluate https://github.com/benruijl/symbolica}"
./bin/build_dependencies.sh "$@"
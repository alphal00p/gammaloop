cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:'aa781d7'}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
./bin/build_dependencies.sh "$@"

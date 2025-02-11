cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/benruijl/symbolica}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b no_namespace https://github.com/alphal00p/symbolica}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
#export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-e534d9f7f8972e22d2a4fb7cd6cb5943373d3bb3}"
./bin/build_dependencies.sh "$@"

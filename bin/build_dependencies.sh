cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/benruijl/symbolica}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b no_namespace https://github.com/alphal00p/symbolica}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
#export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-3f1a5d17abde15df49ca092f311fcd685c1d56fc}"
./bin/build_dependencies.sh "$@"

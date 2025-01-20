cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/benruijl/symbolica}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b nested_symmetrize https://github.com/alphal00p/symbolica}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
#export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-88de4411239548f997dba18cd6a3929d9dc5309c}"
./bin/build_dependencies.sh "$@"

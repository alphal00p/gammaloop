cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}"
#export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-bbdfde9e0df6b6df0699e3ee080c66b3f9b581b8}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
./bin/build_dependencies.sh "$@"

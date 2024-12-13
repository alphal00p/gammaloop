cd ./python/gammaloop
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-b6bc84f640e680d5dd2ca721d0c1ad82d28a2698}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
# export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
./bin/build_dependencies.sh "$@"

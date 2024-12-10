cd ./python/gammaloop
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-371df9de1b1e055bfcb06949c083c037a2fca958}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
# export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
./bin/build_dependencies.sh "$@"

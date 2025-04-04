cd ./python/gammaloop
#export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/benruijl/symbolica}"
# export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b main https://github.com/benruijl/symbolica}"
export CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b symbolica_fork_for_v0_3_3 https://github.com/alphal00p/symbolica}"
export SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-latest}"
#export SYMBOLICA_REVISION_HASH="${SYMBOLICA_REVISION_HASH:-eeed330c824e6ad8f4bc78472864e6d8121b1a12}"
./bin/build_dependencies.sh "$@"

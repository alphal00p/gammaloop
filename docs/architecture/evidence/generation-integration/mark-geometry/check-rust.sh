set -uo pipefail
source target/mark-geometry/environment.sh
cargo nextest run --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --no-fail-fast > target/mark-geometry/geometry-native.log 2>&1
printf '%s\n' "$?" > target/mark-geometry/geometry-native.exit
cargo clippy --locked -p kurvst -p linnest -p clinnet --features linnest/svg --all-targets --profile dev-optim -- --deny warnings > target/mark-geometry/clippy.log 2>&1
printf '%s\n' "$?" > target/mark-geometry/clippy.exit
just fmt --check > target/mark-geometry/format.log 2>&1
printf '%s\n' "$?" > target/mark-geometry/format.exit
cat target/mark-geometry/geometry-native.exit target/mark-geometry/clippy.exit target/mark-geometry/format.exit

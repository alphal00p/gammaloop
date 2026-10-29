set -uo pipefail
source target/mark-geometry-corrected/environment.sh
python3 target/mark-geometry-corrected/native-probe.py > target/mark-geometry-corrected/native-probe.log 2>&1
echo "$?" > target/mark-geometry-corrected/native-probe.exit
cargo nextest run --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --no-fail-fast > target/mark-geometry-corrected/geometry-native.log 2>&1
echo "$?" > target/mark-geometry-corrected/geometry-native.exit
cargo clippy --locked -p kurvst -p linnest -p clinnet --features linnest/svg --all-targets --profile dev-optim -- --deny warnings > target/mark-geometry-corrected/clippy.log 2>&1
echo "$?" > target/mark-geometry-corrected/clippy.exit
just fmt --check > target/mark-geometry-corrected/format.log 2>&1
echo "$?" > target/mark-geometry-corrected/format.exit

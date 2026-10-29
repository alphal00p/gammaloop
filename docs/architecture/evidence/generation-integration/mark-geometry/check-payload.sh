set -uo pipefail
source target/mark-geometry/environment.sh
nix build .#checks.x86_64-linux.linnest-wasm --out-link target/mark-geometry/wasm-check > target/mark-geometry/wasm-check.log 2>&1
echo "$?" > target/mark-geometry/wasm-check.exit
python3 tests/packaging/test_drawing_notices.py -v > target/mark-geometry/packaging.log 2>&1
echo "$?" > target/mark-geometry/packaging.exit
cargo nextest run --locked -p clinnet --cargo-profile dev-optim --no-fail-fast > target/mark-geometry/clinnet.log 2>&1
echo "$?" > target/mark-geometry/clinnet.exit
python3 crates/linnest/layout/native/test_wasm.py -v > target/mark-geometry/runtime.log 2>&1
echo "$?" > target/mark-geometry/runtime.exit
cat target/mark-geometry/wasm-check.exit target/mark-geometry/packaging.exit target/mark-geometry/clinnet.exit target/mark-geometry/runtime.exit

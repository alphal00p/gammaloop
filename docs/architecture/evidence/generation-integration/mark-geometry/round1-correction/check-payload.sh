set -uo pipefail
source target/mark-geometry-corrected/environment.sh
python3 tests/packaging/test_drawing_notices.py -v > target/mark-geometry-corrected/packaging.log 2>&1
echo "$?" > target/mark-geometry-corrected/packaging.exit
cargo nextest run --locked -p clinnet --cargo-profile dev-optim --no-fail-fast > target/mark-geometry-corrected/clinnet.log 2>&1
echo "$?" > target/mark-geometry-corrected/clinnet.exit
python3 crates/linnest/layout/native/test_wasm.py -v > target/mark-geometry-corrected/runtime.log 2>&1
echo "$?" > target/mark-geometry-corrected/runtime.exit
typst compile --root . crates/kurvst/typst/tests/mark-geometry.typ target/mark-geometry-corrected/new-public-fixture.pdf > target/mark-geometry-corrected/new-public-fixture.log 2>&1
echo "$?" > target/mark-geometry-corrected/new-public-fixture.exit
typst query --root . .local/review-mark-geometry-a/boundaries.typ '<boundary-results>' --field value > target/mark-geometry-corrected/new-boundaries-wasm.json 2> target/mark-geometry-corrected/new-boundaries-wasm.log
echo "$?" > target/mark-geometry-corrected/new-boundaries-wasm.exit

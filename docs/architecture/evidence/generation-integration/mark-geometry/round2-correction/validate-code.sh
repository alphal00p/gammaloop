set -u
source docs/architecture/evidence/generation-integration/mark-geometry/round1-correction/environment.sh
out=target/mark-geometry-round3
cargo clippy --locked -p kurvst -p linnest -p clinnet --features linnest/svg --all-targets --profile dev-optim -- --deny warnings > "$out/clippy.log" 2>&1
echo $? > "$out/clippy.exit"
just fmt --check > "$out/format.log" 2>&1
echo $? > "$out/format.exit"
just check-symbolica-feature-isolation > "$out/kernel.log" 2>&1
echo $? > "$out/kernel.exit"
python3 tests/packaging/test_drawing_notices.py -v > "$out/packaging.log" 2>&1
echo $? > "$out/packaging.exit"

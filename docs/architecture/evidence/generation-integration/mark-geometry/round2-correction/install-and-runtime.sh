set -u
source docs/architecture/evidence/generation-integration/mark-geometry/round1-correction/environment.sh
out=target/mark-geometry-round3
bundle=$(readlink -f "$out/clean-result")
# Canonical sync-drawing-assets install/cmp steps, using the already built
# clean commit-addressed canonical package rather than rebuilding a path flake.
install -m 0644 "$bundle/kurvst.wasm" crates/kurvst/typst/kurvst.wasm
install -m 0644 "$bundle/linnest.wasm" crates/linnest/typst/linnest.wasm
cmp "$bundle/kurvst.wasm" crates/kurvst/typst/kurvst.wasm
cmp "$bundle/linnest.wasm" crates/linnest/typst/linnest.wasm
python3 - <<'PY'
from pathlib import Path
p=Path('result');assert p.is_symlink();p.unlink();p.symlink_to((Path('target/mark-geometry-round3/clean-result')).resolve())
PY
typst query --root . .local/review-mark-geometry-round2-a-astra/adversarial-public.typ '<boundary-results>' --field value > "$out/new-public.json" 2> "$out/new-public.log"
echo $? > "$out/new-public.exit"
typst compile --root . crates/kurvst/typst/tests/mark-geometry.typ "$out/new-public-fixture.pdf" > "$out/new-public-fixture.log" 2>&1
echo $? > "$out/new-public-fixture.exit"
python3 crates/linnest/layout/native/test_wasm.py -v > "$out/runtime.log" 2>&1
echo $? > "$out/runtime.exit"
cargo nextest run --locked -p clinnet --cargo-profile dev-optim --no-fail-fast > "$out/clinnet.log" 2>&1
echo $? > "$out/clinnet.exit"
cargo nextest list --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --list-type binaries-only --message-format json > "$out/native-test-binaries.json" 2> "$out/native-test-binaries.log"
cargo nextest list --locked -p clinnet --cargo-profile dev-optim --list-type binaries-only --message-format json > "$out/clinnet-binaries.json" 2> "$out/clinnet-binaries.log"

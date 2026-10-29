set -u
if [ -z "${SYMBOLICA_LICENSE:-}" ]; then export SYMBOLICA_LICENSE="$(cat "${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license")"; fi
export CARGO_BUILD_JOBS=2 NEXTEST_TEST_THREADS=2 RAYON_NUM_THREADS=2 INSTA_UPDATE=no
nix develop /nix/store/b53dxkf559x6hivvnamqs4kzh91d8ar5-nix-shell-env.drv -c bash -c '
set -eu
cargo nextest list --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --message-format json > target/wp3-source38/native-binaries.json 2> target/wp3-source38/native-build.log
cargo nextest run --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --no-fail-fast > target/wp3-source38/native.log 2>&1
cargo test --locked --manifest-path vendor/kurbo/Cargo.toml --lib -- --test-threads=2 > target/wp3-source38/kurbo.log 2>&1
cargo nextest list --locked -p clinnet --cargo-profile dev-optim --message-format json > target/wp3-source38/clinnet-binaries.json 2> target/wp3-source38/clinnet-build.log
cargo nextest run --locked -p clinnet --cargo-profile dev-optim --no-fail-fast > target/wp3-source38/clinnet.log 2>&1
cargo clippy --locked -p kurvst -p linnest -p clinnet --features linnest/svg --all-targets --profile dev-optim -- --deny warnings > target/wp3-source38/clippy.log 2>&1
just fmt --check > target/wp3-source38/format.log 2>&1
just check-symbolica-feature-isolation > target/wp3-source38/kernel.log 2>&1
cargo hakari generate > target/wp3-source38/hakari-generate.log 2>&1
cargo hakari verify > target/wp3-source38/hakari-verify.log 2>&1
just ci-update > target/wp3-source38/ci-update.log 2>&1
python3 tests/packaging/test_drawing_notices.py -v > target/wp3-source38/packaging.log 2>&1
python3 crates/linnest/layout/native/test_wasm.py -v > target/wp3-source38/runtime.log 2>&1
cargo metadata --offline --locked --format-version 1 > target/wp3-source38/root-metadata.json
cargo metadata --offline --locked --manifest-path examples/notebooks/symbolica-host/Cargo.toml --format-version 1 > target/wp3-source38/host-metadata.json
'
rc=$?
printf '%s\n' "$rc" > target/wp3-source38/native.exit
exit "$rc"

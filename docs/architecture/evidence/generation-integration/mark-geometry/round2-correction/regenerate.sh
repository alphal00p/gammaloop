set -u
source docs/architecture/evidence/generation-integration/mark-geometry/round1-correction/environment.sh
out=target/mark-geometry-round3
cargo hakari generate > "$out/hakari-generate.log" 2>&1
echo $? > "$out/hakari-generate.exit"
cargo hakari verify > "$out/hakari-verify.log" 2>&1
echo $? > "$out/hakari-verify.exit"
just ci-update > "$out/ci-update.log" 2>&1
echo $? > "$out/ci-update.exit"

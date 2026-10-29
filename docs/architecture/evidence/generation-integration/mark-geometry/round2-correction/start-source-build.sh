set -euo pipefail
source docs/architecture/evidence/generation-integration/mark-geometry/round1-correction/environment.sh
revision=$(jj log -r @ --no-graph -T commit_id)
flake="git+file:///shared/lcnbr/dev/gammaloop/main?rev=$revision&shallow=1&allRefs=1"
printf '%s\n' "$revision" > target/mark-geometry-round3/build-source-commit
printf '%s\n' "$flake" > target/mark-geometry-round3/clean-flake
jj bookmark create generation-mark-geometry-round3-producer -r "$revision"
nix build "$flake#linnest-wasm" --out-link target/mark-geometry-round3/clean-result > target/mark-geometry-round3/clean-build.log 2>&1
nix derivation show "$flake#linnest-wasm" > target/mark-geometry-round3/clean-derivation.json 2> target/mark-geometry-round3/clean-derivation.log
nix path-info --json "$(readlink -f target/mark-geometry-round3/clean-result)" > target/mark-geometry-round3/clean-package.json
nix build "$flake#checks.x86_64-linux.linnest-wasm" --out-link target/mark-geometry-round3/clean-check > target/mark-geometry-round3/clean-check.log 2>&1
printf '0\n' > target/mark-geometry-round3/clean-build.exit

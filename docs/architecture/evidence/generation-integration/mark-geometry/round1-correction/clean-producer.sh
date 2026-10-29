set -euo pipefail
source target/mark-geometry-corrected/environment.sh
revision=$(cat target/mark-geometry-corrected/versioned-input-commit)
flake="git+file:///shared/lcnbr/dev/gammaloop/main?rev=$revision&shallow=1&allRefs=1"
printf '%s\n' "$flake" > target/mark-geometry-corrected/clean-flake
nix build "$flake#linnest-wasm" --out-link target/mark-geometry-corrected/clean-result > target/mark-geometry-corrected/clean-build.log 2>&1
echo "$?" > target/mark-geometry-corrected/clean-build.exit
nix derivation show "$flake#linnest-wasm" > target/mark-geometry-corrected/clean-derivation.json 2> target/mark-geometry-corrected/clean-derivation.log
nix path-info --json "$(readlink -f target/mark-geometry-corrected/clean-result)" > target/mark-geometry-corrected/clean-package.json 2> target/mark-geometry-corrected/clean-package.log
nix build "$flake#checks.x86_64-linux.linnest-wasm" --out-link target/mark-geometry-corrected/clean-check > target/mark-geometry-corrected/clean-check.log 2>&1
echo "$?" > target/mark-geometry-corrected/clean-check.exit
nix build "$flake#linnest-wasm" --rebuild --keep-failed --no-link > target/mark-geometry-corrected/clean-reproduction.log 2>&1
echo "$?" > target/mark-geometry-corrected/clean-reproduction.exit

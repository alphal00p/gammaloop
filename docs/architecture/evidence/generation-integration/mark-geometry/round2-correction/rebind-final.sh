set -euo pipefail
source target/mark-geometry-round3/environment.sh
out=target/mark-geometry-round3
mkdir -p "$out/intermediate-0bc"
for file in clean-derivation.json clean-package.json clean-build.log clean-check.log clean-reproduction.log clean-reproduction.exit clean-flake build-source-commit; do cp "$out/$file" "$out/intermediate-0bc/$file"; done
nix-store --add-root "$PWD/.local/mark-geometry-round3-intermediate-wasm" --indirect --realise "$(readlink -f "$out/clean-result")" > "$out/intermediate-0bc/retained-store.log"
revision=f35643bae1a67fb6ec0c30151e3624f661cf5f96
flake="git+file:///shared/lcnbr/dev/gammaloop/main?rev=$revision&shallow=1&allRefs=1"
printf '%s\n' "$revision" > "$out/build-source-commit"
printf '%s\n' "$flake" > "$out/clean-flake"
nix build "$flake#linnest-wasm" --out-link "$out/final-result" > "$out/clean-build.log" 2>&1
nix derivation show "$flake#linnest-wasm" > "$out/clean-derivation.json" 2> "$out/clean-derivation.log"
nix path-info --json "$(readlink -f "$out/final-result")" > "$out/clean-package.json"
nix build "$flake#checks.x86_64-linux.linnest-wasm" --out-link "$out/final-check" > "$out/clean-check.log" 2>&1
nix build "$flake#linnest-wasm" --rebuild --keep-failed --no-link > "$out/clean-reproduction.log" 2>&1
printf '0\n' > "$out/clean-reproduction.exit"
cmp "$out/final-result/kurvst.wasm" crates/kurvst/typst/kurvst.wasm
cmp "$out/final-result/linnest.wasm" crates/linnest/typst/linnest.wasm
printf '0\n' > "$out/final-payload-comparison.exit"

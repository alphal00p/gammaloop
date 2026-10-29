set -uo pipefail
source target/mark-geometry/environment.sh
nix derivation show .#linnest-wasm > target/mark-geometry/derivation.json 2> target/mark-geometry/derivation.log
nix path-info --json result > target/mark-geometry/package.json 2> target/mark-geometry/package.log
nix build .#linnest-wasm --rebuild --keep-failed --no-link > target/mark-geometry/wasm-reproduction.log 2>&1
echo "$?" > target/mark-geometry/wasm-reproduction.exit
tail -14 target/mark-geometry/wasm-reproduction.log

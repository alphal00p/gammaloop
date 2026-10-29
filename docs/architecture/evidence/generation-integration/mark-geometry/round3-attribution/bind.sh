set -euo pipefail
source target/mark-geometry-attribution/environment.sh
revision=$(jj log -r @ --no-graph -T commit_id)
flake="git+file:///shared/lcnbr/dev/gammaloop/main?rev=$revision&shallow=1&allRefs=1"
nix eval --json "$flake#packages.x86_64-linux.linnest-wasm" --apply 'p: { source = p.src.outPath; derivation = p.drvPath; }' > target/mark-geometry-attribution/runtime-source.json 2> target/mark-geometry-attribution/runtime-source.log
nix hash path /nix/store/88qva12y43c0fig08mslcrmjbb3ydxa9-linnest-wasm-0.1.0 > target/mark-geometry-attribution/nar-hash.txt

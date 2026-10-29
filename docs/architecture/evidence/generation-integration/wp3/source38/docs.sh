set -eu
nix develop /nix/store/b53dxkf559x6hivvnamqs4kzh91d8ar5-nix-shell-env.drv -c bash -c '
set -eu
python3 target/wp3-source38/audit.py > target/wp3-source38/audit.log
for note in generation-mark-geometry generation-integration-plan; do
 typst compile --features html --root . --creation-timestamp 0 docs/architecture/$note.typ target/wp3-source38/$note.pdf > target/wp3-source38/$note-pdf.log 2>&1
 printf "#document(\"index.html\")[#include \"/docs/architecture/%s.typ\"]\n" "$note" > target/wp3-source38/$note-wrapper.typ
 typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 target/wp3-source38/$note-wrapper.typ target/wp3-source38/$note-html > target/wp3-source38/$note-html.log 2>&1
done
maturin sdist --manifest-path crates/gammaloop-api/Cargo.toml --out target/wp3-source38/sdist > target/wp3-source38/sdist.log 2>&1
'
printf '0\n' > target/wp3-source38/docs.exit

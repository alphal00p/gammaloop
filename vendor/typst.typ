= Typst document compiler dependencies

The documentation builder's optional `watch` feature uses upstream Typst 0.15.1.
Its internal component family remains pinned in `Cargo.lock`, matching the
Nix-provided CLI used for release output parity. No Typst compiler component
requires a local patch.

The removed `typst-library` and `typst-svg` patches only introduced optional
feature gates for the former notebook compiler: plugin loading, syntax
highlighting, bibliography, PDF images, raster formats and PDF-to-SVG conversion.
Comparison with the published 0.15.1 sources showed those gates and their
disabled-feature substitutes were the only implementation changes. The
builder enabled every gate (`typst-library/full` and `typst-svg/pdf`), restoring
upstream behavior. Its direct dependencies on these components existed only
to select those patched features; the upstream `typst` dependency owns them.

Notebook math and graph labels render directly as MathML with native Rust
geometry. The Python runtime does not link this Rust compiler; explicit vector
exports use the optional Python `typst` package. Graph layout and line drawing
never call a Typst plugin.

== Fonts

The documentation builder uses Typst Kit's embedded and system fonts.
The independent browser STIX font and its license are owned by
`linnet-render-py` and included in `gammaloop-render-data`. They are discovered
as resources rather than bundled with the notebook host extension.

== Updating and validation

Keep the document compiler's component versions aligned when upgrading Typst,
regenerate the Cargo/Nix graph, and run the documentation builder's checks and
tests, including its `watch` feature. Runtime browser checks separately exercise
MathML, measurement receipts and foreignObject SVG labels. Check resolved target
dependencies to ensure document compiler crates remain outside runtime builds.

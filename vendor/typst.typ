= Typst document compiler features

`typst-library` and `typst-svg` are vendored from the published crates.io
0.15.1 sources. Their Apache-2.0 license is retained in each directory.
The registry metadata, lockfiles, and workspace-dependent `Cargo.toml.orig`
are omitted; normalized Cargo manifests remain buildable outside upstream's
workspace. Modified files identify this patch in their first line.

Notebook math and graph labels now render directly as MathML. The Python
runtime does not link this Rust compiler; explicit vector exports use the
optional Python `typst` package. The patches below remain for document-building
tools and the standalone compiler crate. Graph layout and line drawing run
directly in Rust and never call a Typst plugin. The compiler library does
not expose upstream feature flags for several large document-only components,
so this patch introduces opt-in features with empty defaults:

- `typst-library/plugins`: WebAssembly plugin loader and Wasmi interpreter.
- `typst-library/syntax-highlighting`: Syntect and its Two Face syntax archive.
  Without it, raw text retains its whitespace and layout but has no coloring;
  explicitly loading a theme or syntax returns a diagnostic.
- `typst-library/bibliography`: Hayagriva, CSL styles, BibLaTeX, and collation.
  Without it, citation and bibliography constructors are absent. Shared layout
  records remain available to upstream HTML/layout crates, but no citation
  processor or rendered works can be constructed.
- `typst-library/pdf`: PDF image parsing. Without it, loading a PDF image
  returns a diagnostic; the shared image type is uninhabited.
- `typst-library/raster-formats`: JPEG, GIF, and WebP decoders. PNG/raw pixels
  remain for SVG and font image payloads.
- `typst-svg/pdf`: PDF-to-SVG conversion, including Hayro and its standard fonts.
  Enable this together with library PDF support when using the SVG exporter.
- `typst-library/full`: all library features above.

`typst-renderer` uses this lean compiler by default. Its `full` feature enables
all the above plus PDF/PNG export and Typst Kit's embedded/system fonts.
The persistent documentation builder selects full library and SVG support.
Workspaces compiling these document tools must carry both Cargo patches. The
isolated notebook host no longer needs either patch. Cargo features are additive,
so document-tool builds and runtime dependency checks remain separate.

The math runtime retains the rest of Typst's shared evaluator and layout code;
this is not a separate implementation of Typst math. Do not remove shaping,
Unicode segmentation, SVG parsing, or MathML support merely because graph
geometry uses direct SVG.

== Fonts

`crates/typst-renderer/fonts` holds seven unmodified fonts from
`typst-assets` 0.15.1: Libertinus Serif regular/bold/italic/bold-italic,
New Computer Modern Math Book/bold, and DejaVu Sans Mono regular.
These cover document math, text styles, and raw text. Typst equations request
weight 450, so the Book face preserves the upstream default math appearance. Keeping complete math
fonts preserves arbitrary user labels; no glyph subsetting is applied.
The selected files total 4,227,456 bytes, replacing the 17-font bundle.
Their upstream `NOTICE` is preserved alongside them and copied into the
document asset bundle. Update that copy when changing the fonts. The browser
runtime's independent STIX font is owned by `linnet-render-py`.

== Updating and validation

Compare these directories against the published 0.15.1 crate sources to isolate
the feature patch. On upgrade, reapply the feature boundaries, update the root
workspace lockfile, regenerate the Cargo/Nix graph, and run both lean and full
`typst-renderer` tests. The full test exercises syntax highlighting, citations,
PDF/PNG export, and PDF-to-SVG import; lean tests exercise math, text, SVG import,
MathML, and diagnostics for omitted features. Runtime browser checks instead
exercise direct MathML, measurement receipts, and foreignObject SVG labels.
Inspect resolved target dependencies to ensure these document compiler crates
have not returned through runtime feature unification.

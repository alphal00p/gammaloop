= Typst notebook compiler features

`typst-library` and `typst-svg` are vendored from the published crates.io
0.15.1 sources. Their Apache-2.0 license is retained in each directory.
The registry metadata, lockfiles, and workspace-dependent `Cargo.toml.orig`
are omitted; normalized Cargo manifests remain buildable outside upstream's
workspace. Modified files identify this patch in their first line.

`typst-renderer` names both directories as path dependencies as well as the
workspace patches, so Maturin carries them into the `linnet-py` source
distribution. Each directory's added `rustfmt.toml` keeps `cargo fmt --all`
from reformatting the vendored source. The other `vendor/` crates follow the
same convention.

The notebook needs Typst's language evaluator, text shaping, math layout,
SVG exporter/importer, and HTML/MathML exporter. Native graph drawing lays out
and draws lines directly in Rust without a Typst plugin; Typst figures and
templates load the CeTZ, Linnest, Kurvst, and MiTeX packages, which call
plugins. The compiler library does not expose upstream feature flags for
several large document-only components, so this patch introduces opt-in
features with empty defaults:

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
all the above plus PDF/PNG export; either way, the fonts are the embedded set
below. `linnet-py` selects `full`, so its standalone wheel and the notebook host
that links it carry every feature; the persistent documentation builder selects
full library and SVG support. Consumer workspaces must carry both Cargo
patches, as the repository and isolated notebook host do. Cargo features are
additive: build the notebook from its isolated host manifest so that other
workspace members cannot add features to its dependencies.

The math runtime retains the rest of Typst's shared evaluator and layout code;
this is not a separate implementation of Typst math. Do not remove shaping,
Unicode segmentation, SVG parsing, or MathML support merely because graph
geometry uses direct SVG.

== Fonts

`crates/typst-renderer/fonts` holds eleven unmodified fonts, all under the SIL
Open Font License or the Bitstream Vera license: Libertinus Serif regular,
semibold, and bold, each upright and italic, and DejaVu Sans Mono regular, bold,
oblique, and bold oblique from `typst-assets` 0.15.1, and STIX Two Math 2.13
b171 from the STIX Fonts project. They cover text styles, math, and raw text and
total 3,874,184 bytes. They are the only fonts. Do not enable Typst Kit's
`embedded-fonts` feature in a distributed build: its set adds the New Computer
Modern fonts, of which NewCM10-Regular is under the GPL (with a font exception)
and the others under the GUST Font License. Keeping complete fonts preserves
arbitrary user labels; no glyph subsetting is applied.

Typst's default math font, New Computer Modern Math, is not bundled.
`typst-renderer` therefore adds a library show-set rule that typesets equations
in STIX Two Math; without it, equations fall back to a text font and lose their
math layout. STIX Two Math has a single weight, so bold text does not embolden
equations; `bold(..)` in math still selects bold glyphs.

`fonts/NOTICE` holds the typst-assets notice, without the sections for fonts
that are not distributed, and the STIX Two Math notice. It is copied into the
license directories of the notebook wheel and `linnet-py`. Update those copies
when changing the fonts.

== Updating and validation

Compare these directories against the published 0.15.1 crate sources to isolate
the feature patch. On upgrade, reapply the feature boundaries, update both
workspace lockfiles, regenerate the Cargo/Nix graph, and run the
`typst-renderer` tests without features and with `full`. The full test
exercises plugins, syntax highlighting, citations, PDF/PNG export, and
PDF-to-SVG import; the lean tests exercise diagnostics for omitted features;
both exercise the embedded font variants, math, text, SVG import, and MathML.
Then build the isolated Pyodide wheel and run `installed_offline_rendering.py`
and the tensor network display checks in a browser. Inspect resolved target
dependencies to ensure optional components have not returned through feature
unification.

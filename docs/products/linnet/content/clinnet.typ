#import "../../shared.typ": callout, boundary, source-link

#let clinnet = [
= Render DOT figures with Clinnet

Clinnet is Linnet's command-line companion for documentation and diagram workflows. The Cargo
package is named `clinnet`, while the installed executable is named `linnet`. It recursively
finds DOT files, renders one PDF per graph through Typst, and assembles those PDFs into a grid.
Use the Rust library for graph algorithms, the `linnet-py` `Graph` API for typed programmatic
drawing, and Clinnet when the input is already DOT and the desired result is a browsable set of
figures.

== Install and render a directory

Clinnet requires Typst 0.15.0 or newer. Install both commands, then confirm the executable's
current flags rather than relying on an older pre-subcommand invocation:

// docs-example: syntax clinnet-render-dot
```sh
cargo install clinnet
cargo install typst-cli --version 0.15.0 --locked
linnet --help
linnet draw examples
```

The final command scans `examples/` recursively. With default paths it creates:

- `build/figs/`, mirroring the input directory and containing one PDF per DOT file;
- `build/grid.pdf`, the combined document;
- `build/templates/`, containing the default figure, grid, layout, Linnest, and Kurvst assets;
- `build/.cache/figures.json`, the per-figure content cache; and
- `build/.cache/run-metadata.json`, the plan used by targeted redraw commands.

#callout("What a second run reuses", [
  A figure is reused only when its DOT source, selected figure template, explicit `--style`
  files, `--input` values, and the bundled `.typ`/`.wasm` assets in `build/templates/` have the
  same content hash. Clinnet does not recursively discover arbitrary imports beside a custom
  template; pass those files with `--style` when they must invalidate the cache. The grid is
  still regenerated from the current figure index. Removing a DOT file removes its stale figure
  output during the next full `draw`; deleting the figure cache forces every remaining figure to
  rebuild.
])

== Configure the V1 templates

Use `--figure-template` for each individual figure and `--grid-template` for the combined PDF.
Repeat `--style` for extra files whose contents must invalidate the figure cache. Repeat
`--input key=value` for the CLI's recognized typed configuration overrides:

// docs-example: syntax
```sh
linnet --build-dir build draw graphs \
  --figure-template templates/figure.typ \
  --grid-template templates/grid.typ \
  --style templates/colors.typ \
  --input steps=250 \
  --input theme=dark \
  --input accent=teal \
  --input columns=4 \
  --output-path build/review.pdf
```

Clinnet parses values before constructing the native Typst configuration. Figure `steps` and
`seed` become generic layout options; every other accepted figure input is staged as typed
`config.options` for the selected template to interpret. In the example, `theme` and `accent`
are application-defined settings rather than Clinnet modes. Grid overrides cover its rows,
columns, page geometry, alignment, labels, and page numbers. Inputs are never published as
`sys.inputs`.

Each figure render stages the DOT topology and an ephemeral entrypoint. That entrypoint imports
the selected template, constructs `config`, adds `data-path`, and calls the mandatory
`render(config)` export:

```typst
#import "layout.typ": layout

#let render(config) = {
  set page(width: auto, height: auto)
  context layout(config)
}
```

The Typst command receives only the entrypoint, output, and root paths, so large configurations
do not produce large process argument lists. `--output-path` selects the final grid PDF; it is not
a figure configuration field. Plain strings stay escaped data, and custom templates do not have a
legacy top-level compilation path. The embedded defaults include the canonical Linnest and Kurvst
package trees, so a custom template can use the same graph-layout and curve primitives. Continue
to the
#link("guides/linnest/")[Linnest Typst API] for graph construction, layout, inspection, and
drawing rather than copying its wrapper surface into a Clinnet template guide.

#callout("Generic and GammaLoop templates", [
  Clinnet's bundled template is generic: it attaches typed element records, applies final labels
  and styles before ordered layout passes, and draws the result. GammaLoop writes a separate V1
  bundle under `drawings/templates/`, including its generated `edge-style.typ`. Select that
  bundle's `figure.typ` when particle colors and decorations, momentum annotations, amplitude
  ordering, or cut-matched cross-section labels are required. Both templates implement the same
  `render(config)` contract, while only GammaLoop's template assigns meaning to its domain
  options.
])

`--style` records an imported file as a cache dependency; it does not inject or execute that file.
The selected template must import it normally through Typst's module system. Clinnet does not
recursively discover every custom import, so repeat `--style` for imports whose changes must
invalidate cached figures.

#boundary("Templates must share a Typst project root", [
  Clinnet chooses the narrowest existing directory that contains the template, staged topology,
  and explicit renderer dependencies. A Typst import outside that tree, a missing style file, or
  a template moved without its dependencies fails before a PDF is produced. Keep templates and
  imported assets beneath one project directory and read the full Typst error printed by
  Clinnet.
])

== Make a targeted rebuild

Run a complete `draw` once so that the build metadata records the input root, templates, output
paths, styles, and figures. You can then rebuild one known DOT figure without rescanning the
directory, or rebuild only the grid:

// docs-example: syntax
```sh
linnet --build-dir build redraw-figure graphs/box.dot --input theme=light
linnet --build-dir build redraw-grid --output-path build/review.pdf
```

`redraw-figure` starts with the inputs recorded by the full draw and applies new inputs as
per-key overrides. `redraw-grid` reuses the cached figure list; if any `--input` is supplied to
that command, the supplied list replaces the grid's recorded input list, so repeat every value
the grid template still needs. A target not present in the preceding plan must be added with a
new full `draw`.

Use `--no-parallel` for deterministic sequential diagnosis, or `--jobs N` to bound parallel
figure compilation. These options affect scheduling, not cache identity or graph semantics.

== Source and release boundaries

The implementation and exact Clap definitions are in
#source-link("crates/clinnet/src/main.rs", label: "the Clinnet command source"). Clinnet has its
own package version and #link("version-history/clinnet/")[version history], independent of the
Rust `linnet` library. Record the Clinnet and Typst versions when reporting a rendering or cache
problem.
]

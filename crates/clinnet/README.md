# clinnet

`clinnet` is the command-line companion to the [`linnet`](https://crates.io/crates/linnet) graph library.
It scans a directory for `.dot` files, incrementally renders each figure through
[Typst](https://typst.app), and assembles the resulting PDFs into a single grid document.

## Installation

```bash
cargo install clinnet
cargo install typst-cli --version 0.15.0 --locked
```

This installs the `linnet` binary (the package is named `clinnet` to avoid a clash with the library
crate). Linnet requires Typst 0.15.0 or newer and shells out to the `typst` CLI for every
compilation step.

## Quick start

```bash
# Render every .dot file underneath ./examples into build/grid.pdf
linnet draw examples
```

On first run, clinnet writes its default Typst templates and the Linnest/Kurvst package trees into
`build/templates/`. You can provide your own templates via `--figure-template` and
`--grid-template`. Pass `--style` multiple times to add extra files whose contents should
invalidate the incremental figure cache when they change.

You can pass custom variables to Typst templates using `--input key=value`, which makes the
variables available through `sys.inputs` in your templates. This works exactly like the Typst CLI's
`--input` option.

## Typical workflow

1. Generate or edit `.dot` graphs in your project directory.
2. Run `linnet draw <graphs-root>` to refresh the figure PDFs and the combined grid.
3. Open `build/grid.pdf` (or the path supplied via `--output-path`) to browse the results.

The cache lives in `build/.cache/figures.json`; deleting it forces a full rebuild. See
`linnet --help` for the complete flag list.

## Passing variables to templates

You can pass custom variables to your Typst templates:

```bash
# Pass variables that will be available in templates via sys.inputs
linnet draw examples --input theme=dark --input scale=1.5

# Multiple inputs can be specified
linnet draw examples --input author="John Doe" --input version=1.0
```

These variables are available in both figure and grid templates through Typst's `sys.inputs`
dictionary.

## Typst wrapper API

The generated templates include the canonical Linnest Typst package tree and its
WASM plugin under `crates/linnest/typst`. Import the package source relative to
the generated template directory:

```typst
#import "crates/linnest/typst/src/lib.typ": draw, graph, layout, subgraph

#import graph: build, edge, node, sink, source

#let g = build({
  node(<a>, label: [A])
  node(<b>, label: [B])
  edge(
    source(<a>, compass: "e"),
    <a-b>,
    sink(<b>, compass: "w"),
    label: [a-b],
  )
}, name: "example")
#let g = layout(g)
#let info = graph.info(g)
#let nodes = graph.nodes(g)
#let north = subgraph.compass(g, "n")
#let north-edges = graph.edges(g, subgraph: north)
#let dot = graph.dot(g)
#draw(g)
```

The wrapper exports:

- `draw(graph, ...)`: draws a laid-out graph object through CeTZ.
- `layout(graph, seed: 2, steps: 5, ...)`: lays out one graph object.
- `graph`: namespace for graph construction, DOT parsing/printing, inspection, joins, and graph algorithms.
- `subgraph`: namespace for subgraph object construction and inspection.

See `crates/linnest/typst/docs/manual.pdf` for the full Tidy-generated API reference, and `crates/linnest/examples/typst-api.typ` for a compile-checked example:

```bash
typst compile --root . crates/linnest/examples/typst-api.typ /tmp/linnest-typst-api.pdf
```

## Partial rebuilds

After a full run, you can re-render individual figures or the final grid without rescanning DOT
files:

```bash
# Rebuild a single figure that was part of the previous run
linnet --build-dir build redraw-figure path/to/file.dot

# Rebuild the grid PDF based on the cached figure metadata
linnet --build-dir build redraw-grid

# Partial rebuilds also support input variables
linnet --build-dir build redraw-grid --input theme=light
```

`redraw-figure` starts with the inputs cached by the preceding `draw` and applies any supplied
`--input` values as per-key overrides.

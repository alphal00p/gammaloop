# clinnet

`clinnet` is the command-line companion to the [`linnet`](https://crates.io/crates/linnet) graph library.
It scans a directory for `.dot` files, incrementally renders each figure through
[Typst](https://typst.app), and assembles the resulting PDFs into a single grid document.

## Installation

```bash
cargo install clinnet
```

This installs the `linnet` binary (the package is named `clinnet` to avoid a clash with the library
crate). Make sure the `typst` CLI is available on your `$PATH`, as the tool shells out to it for
every compilation step.

## Quick start

```bash
# Render every .dot file underneath ./examples into build/grid.pdf
linnet examples
```

On first run, clinnet writes a small set of default Typst templates and a compiled `linnest.wasm`
plugin into `build/templates/`. You can provide your own templates via `--figure-template`,
`--grid-template`, or `--layout-template` flags. Pass `--style` multiple times to add extra files
whose contents should invalidate the incremental build cache when they change.

You can pass custom variables to Typst templates using `--input key=value`, which makes the
variables available through `sys.inputs` in your templates. This works exactly like the Typst CLI's
`--input` option.

## Typical workflow

1. Generate or edit `.dot` graphs in your project directory.
2. Run `linnet <graphs-root>` to refresh the figure PDFs and the combined grid.
3. Open `build/grid.pdf` (or the path supplied via `--grid-output`) to browse the results.

The cache lives in `build/.cache/figures.json`; deleting it forces a full rebuild. See
`linnet --help` for the complete flag list.

## Passing variables to templates

You can pass custom variables to your Typst templates:

```bash
# Pass variables that will be available in templates via sys.inputs
linnet examples --input theme=dark --input scale=1.5

# Multiple inputs can be specified
linnet examples --input author="John Doe" --input version=1.0 --input debug=true
```

These variables are available in both figure and grid templates through Typst's `sys.inputs`
dictionary.

## Typst wrapper API

The generated templates include `linnest.typ`, a small Typst wrapper around `linnest.wasm`.
Import it from a directory that also contains the WASM plugin:

```typst
#import "linnest.typ": (
  builder,
  builder-add-edge,
  builder-add-node,
  finish-builder,
  graph-edges,
  graph-info,
  graph-nodes,
  layout-graph,
  subgraph-from-compass,
)

#let b0 = builder(spec: (
  name: "example",
))
#let a = builder-add-node(b0, name: "a")
#let b = builder-add-node(a.builder, name: "b")
#let b2 = builder-add-edge(b.builder, source: (node: a.node), sink: (node: b.node))
#let graph = finish-builder(b2)
#let graph = layout-graph(graph)
#let info = graph-info(graph)
#let nodes = graph-nodes(graph)
#let north = subgraph-from-compass(graph, "n")
#let north-edges = graph-edges(graph, subgraph: north)
```

The wrapper exports:

- `default-layout-config`: the default layout settings used by `layout.typ`.
- `builder(spec: (:))`: creates an archived graph builder byte array.
- `builder-add-node(builder, name: none, index: none, statements: (:))`: adds a node and returns `(builder: bytes, node: index)`.
- `builder-add-edge(builder, source: none, sink: none, orientation: "default", flow: none, id: none, statements: (:))`: adds a paired or external edge and returns the next builder.
- `finish-builder(builder)`: turns an archived builder into an archived graph.
- `build-graph(spec)`: convenience sugar that builds one archived graph byte array from a Typst dictionary.
- `parse-graphs(input)`: parses one or more DOT digraphs and returns archived graph byte arrays.
- `layout-graph(graph, config: default-layout-config)`: lays out one parsed graph.
- `layout-graphs(input, config: default-layout-config)`: parses and lays out all DOT graphs in `input`.
- `graph-info(graph)`: returns graph metadata and default DOT statements.
- `graph-nodes(graph, subgraph: none)`: returns node records, optionally filtered by an archived subgraph.
- `graph-edges(graph, subgraph: none)`: returns edge records, optionally filtered by an archived subgraph.
- `subgraph-from-label(graph, label)`, `subgraph-from-bits(graph, bits)`, and `subgraph-from-compass(graph, compass)`: construct archived subgraph byte arrays.
- `subgraph-label(subgraph)`, `subgraph-hedges(subgraph)`, and `subgraph-contains-hedge(subgraph, hedge)`: inspect archived subgraphs.
- `graph-cycle-basis(graph)`: returns archived subgraphs for the cycle basis.
- `graph-spanning-forests(graph)`: returns archived subgraphs for spanning forests.
- `join-graphs(left, right, key: "...")`: joins dangling edges whose half-edge data field matches on `key`; supported keys are `statement`, `port_label`, `compass`, and `id`.

See `crates/linnest/examples/typst-api.typ` for a compile-checked example:

```bash
typst compile --root . crates/linnest/examples/typst-api.typ /tmp/linnest-typst-api.pdf
```

## Partial rebuilds

After a full run, you can re-render individual figures or the final grid without rescanning DOT
files:

```bash
# Rebuild a single figure that was part of the previous run
linnet --build-dir build --rebuild-figure path/to/file.dot

# Rebuild the grid PDF based on the cached figure metadata
linnet --build-dir build --rebuild-grid

# Partial rebuilds also support input variables
linnet --build-dir build --input theme=light --rebuild-grid
```

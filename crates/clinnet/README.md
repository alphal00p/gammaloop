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
#import "linnest.typ": parse-graphs, layout-graph, graph-info, graph-nodes, graph-edges

#let graphs = parse-graphs(read("graph.dot"))
#let graph = layout-graph(graphs.first())
#let info = graph-info(graph)
#let nodes = graph-nodes(graph)
#let edges = graph-edges(graph)
```

The wrapper exports:

- `default-layout-config`: the default layout settings used by `layout.typ`.
- `parse-graphs(input)`: parses one or more DOT digraphs and returns archived graph byte arrays.
- `layout-graph(graph, config: default-layout-config)`: lays out one parsed graph.
- `layout-graphs(input, config: default-layout-config)`: parses and lays out all DOT graphs in `input`.
- `graph-info(graph)`: returns graph metadata and default DOT statements.
- `graph-nodes(graph, subgraph: none)`: returns node records, optionally filtered by a subgraph label or bool array.
- `graph-edges(graph, subgraph: none)`: returns edge records, optionally filtered by a subgraph label or bool array.
- `graph-subgraph(graph, subgraph)`: normalizes a subgraph spec into the base62 label used by the Rust API.
- `graph-compass-subgraph(graph, compass)`: builds a base62 subgraph label from a compass point such as `"n"` or `"s"`.

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

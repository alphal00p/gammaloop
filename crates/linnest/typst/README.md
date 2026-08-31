# Linnest

Linnest is the Typst wrapper for the `linnest.wasm` graph-layout plugin. It requires Typst 0.15.0
or newer.

`src/render/figure.typ` owns the domain-neutral `render(config)` implementation used by native
graph-spec consumers. DOT-oriented tools can pass a parsed graph to the shared
`layout-graph(config, graph)` adapter instead of duplicating the styling and layout pipeline.

The canonical source manual is [`docs/manual.typ`](docs/manual.typ). Its workflow is published
as the [Linnest Typst guide](https://alphal00p.github.io/gammaloop/products/linnet/latest/guides/linnest/),
with the supported symbol reference split across the
[Linnest Typst API reference](https://alphal00p.github.io/gammaloop/products/linnet/latest/reference/typst/).
Compile-checked examples live under [`examples/`](examples/).

Contributor policy is in [`CONTRIBUTING.typ`](../../../CONTRIBUTING.typ).

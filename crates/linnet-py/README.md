# Linnet Community bindings

This Rust binding library supplies `symbolica.community.graph` in the complete
Symbolica Community wheel. It is not a standalone Python distribution.

See [the Python guide](../../docs/products/linnet/content/quickstart-python.typ)
and [the API guide](../../docs/products/linnet/content/api.typ) for builders,
topology, graph algorithms, DOT, layout streaming, and rendering.

Run the Python suites in `tests/` against an installed Community wheel. The
repository notebook host also embeds these bindings with HepKit and Tensor.

`render()` returns a shared `DiagramView`. Its default `snapshot()`, `to_svg()`
and `to_html()` produce static output without launching a browser. Install
`anywidget` for live notebook measurement, or explicitly call
`measured(browser="chromium")` with Playwright. `measured()` requires completed
measurements; `outlined()` uses optional typst-py for vector labels and PDF/PNG.
Fonts are discovered through `symbolica.data` and `symbolica-fonts`.
See the [migration notes](../../docs/products/feynkit/content/community-host.typ)
for the MIME display contract and snapshot defaults.

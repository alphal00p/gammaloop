# Linnet Community bindings

This Rust binding library supplies `symbolica.community.graph` in the complete
Symbolica Community wheel. PyO3 fixes each class's module at compile time: the
`community` feature compiles them as `symbolica.community.graph`, and builds
without it as the standalone module `linnet_py`. No standalone wheel is built
yet; `linnet_py.pyi` describes its surface, whose graph keys exclude Symbolica
expressions.

See [the Python guide](../../docs/products/linnet/content/quickstart-python.typ)
and [the API guide](../../docs/products/linnet/content/api.typ) for builders,
topology, graph algorithms, DOT, layout streaming, and rendering.

Run the Python suites in `tests/` against an installed Community wheel. The
repository notebook host also embeds these bindings with HepKit and Tensor.

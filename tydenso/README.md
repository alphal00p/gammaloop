# Tydenso

Tydenso is GammaLoop's Typst interface to Spenso and Idenso. Its nested Cargo
workspace intentionally keeps the WebAssembly dependency graph separate from
GammaLoop's native Python and command-line builds.

The complete Typst package lives in [`typst`](typst), including its manual,
examples, tests, compressed engine, and small inflater plugin. Build and check
it from the repository root with:

```sh
just tydenso::build
TYMBOLICA_CHECKOUT=/path/to/symbolica-typst-plugin just tydenso::check
just tydenso::manual
```

The interop tests require `TYMBOLICA_CHECKOUT` to point to
[`symbolica-dev/symbolica-typst-plugin`](https://github.com/symbolica-dev/symbolica-typst-plugin) at Git HEAD
`cf7b9fe59ba3fc7fe9ea70b875f5660ee2addd7a`, the exact Symbolica Typst plugin revision pinned
by the nested Rust workspace. The check uses Nix to rebuild its combined algebra and integration engine
in a temporary copy. GammaLoop, Tydenso and the Typst plugin must all resolve
the same Symbolica Git revision from the root workspace's patch, so their
Atom payload formats agree; Numerica and Graphica remain at crates.io `3.0.0`.
It preserves the checkout and reuses compiled dependencies in
`target/tymbolica` at the repository root. The Rust payload dependency comes
from the same pinned Git revision.

See the [Tydenso manual source](typst/manual.typ) for the public Typst interface.

Tydenso is licensed under the [MIT license](LICENSE) carried in this directory.

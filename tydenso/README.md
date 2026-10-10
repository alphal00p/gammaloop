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
`b96c1f3ebd1c2563bf91079a0dd6dc5ec82191c0`, its Symbolica 3.0.1 release update.
The check uses Nix to rebuild that engine in a temporary copy. GammaLoop,
Tydenso and the engine resolve Symbolica and Numerica `3.0.1`, and Graphica
`3.0.0`, from crates.io. The separate Rust payload library remains pinned to
`7f869adc14dcf2758bab16aaa25e06eb10ffbfad`; it uses the same envelope implementation
and backend format. The check preserves the checkout and reuses compiled
dependencies in `target/tymbolica` at the repository root.

Tydenso is licensed under the [MIT license](LICENSE) carried in this directory.

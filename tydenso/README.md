# Tydenso

Tydenso is GammaLoop's Typst interface to Spenso and Idenso. Its nested Cargo
workspace intentionally keeps the WebAssembly dependency graph separate from
GammaLoop's native Python and command-line builds.

The complete Typst package lives in [`typst`](typst), including its manual,
examples, tests, compressed engine, and small inflater plugin. Build and check
it from the repository root with:

```sh
just tydenso::build
TYMBOLICA_CHECKOUT=/path/to/tymbolica just tydenso::check
just tydenso::manual
```

The interop tests require `TYMBOLICA_CHECKOUT` to have Git HEAD
`32ae00483b4845161401afac8be755cc336cb2bb`, the exact Tymbolica revision pinned
by the nested Rust workspace. The check uses Nix to rebuild Tymbolica and Rubi
in a temporary copy with the same Symbolica revision as GammaLoop and Tydenso.
It preserves the checkout and reuses compiled dependencies in
`target/tymbolica` at the repository root. The Rust payload dependency comes
from the same pinned Git revision.

Tydenso is licensed under the [MIT license](LICENSE) carried in this directory.

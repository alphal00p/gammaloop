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
`e85233dd3d99e5fb458913d0d49f99d83466025b`, the exact Tymbolica revision pinned
by the nested Rust workspace. The checkout only supplies Tymbolica's Typst
package files to those tests; the Rust payload dependency itself comes from the
same pinned Git revision.

Tydenso is licensed under the [MIT license](LICENSE) carried in this directory.

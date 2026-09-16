#import "../../shared.typ": source-link

#let community-host = [
= Integrating the Symbolica community host

`feynkit-py` is an `rlib` implementing `SymbolicaCommunityModule`. The combined
#link("https://github.com/symbolica-dev/symbolica-community")[Symbolica community host] owns the
single `symbolica.core` extension. Link FeynKit into that host so all community modules exchange
the same Symbolica expression type and global state.

== Register and package the module

Add `feynkit-py` from the checkout revision that the host is testing. During local integration a
path dependency is explicit:

// docs-example: syntax
```toml
[dependencies]
feynkit-py = { path = "../gammaloop/crates/feynkit-py" }
```

For a release, replace the path with the repository URL and a tested Git `rev`. Keep Symbolica,
PyO3, and the other community modules compatible with this checkout. The `ufo` feature is enabled
by default in the Python adapter; disable default features if the host intentionally omits UFO
loading. The pure Rust facade has different defaults.

In the host's existing core module, after `create_symbolica_module(m)?`, use its registration
macro:

// docs-example: syntax
```rust
register_module!(m, feynkit_py::FeynkitModule);
```

The host creates `symbolica.community.feynkit_native` and supplies `initialize_module()`.
Copy the package containing
#source-link("crates/feynkit-py/python/symbolica/community/feynkit/__init__.py", label: "the Python wrapper")
into the host's `python/symbolica/community` tree. Its wrapper imports the native module and
calls the initializer. FeynKit itself must not declare another PyO3 extension entry point.
Include `typst>=0.15,<0.16` in the host's display dependencies for automatic notebook figures;
include `ufo-model-loader` when offering raw UFO import.

== Export the documented Python surface

Forward `feynkit-py/python_stubgen` from the host's stub-generation feature, alongside its other
community modules. The documentation exporter updates both the package stub and the reference
input from the same native module:

// docs-example: syntax
```sh
cargo run --locked -p alphal00p-docs-python-exporter --features feynkit -- \
  feynkit-community docs/api/python/feynkit-community.pyi
cargo run --locked -p alphal00p-docs-python-exporter --features feynkit -- \
  feynkit-community docs/api/python/feynkit-community.pyi --check
```

The public stub is generated from native signatures and docstrings. Its documentation audit
requires examples and parameter descriptions. The product's
#link("reference/python/feynkit-community/")[Python API pages] use that same registered surface;
refresh and check the generated inputs when the API changes. The check also validates runtime
exports and the syntax of documented Python examples.

Build the product documentation or keep its local preview running with:

// docs-example: syntax
```sh
just docs-site feynkit
just docs-watch feynkit 8117
```

The preview is served at `http://127.0.0.1:8117`. These commands use the registered Rust and
Python components, manual pages, and example catalog.

== Build the browser showcases

The #link("guides/showcases/")[showcase gallery] runs in Marimo's browser Python runtime.
The repository's combined host packages Symbolica, FeynKit, Spenso, and Idenso into one
WebAssembly wheel. Its `wasm` feature selects portable numeric backends; its default
`native` feature retains the desktop backends. Browser Symbolica runs without a license key.

The wheel command uses the checkout's pinned Emscripten Rust toolchain and provisions the
matching Pyodide build environment through cibuildwheel. Then export the executable
cells and start the existing documentation watcher:

// docs-example: syntax
```sh
just notebook-wheel
just notebook-ufo-wheel
just docs-notebooks /path/to/symbolica-wasm.whl feynkit docs/generated/notebooks target/notebook-ufo-wheel/ufo_model_loader-0.1.8-py3-none-any.whl
just docs-watch feynkit 8117
```

Pass the actual Symbolica wheel filename produced by the build. The UFO loader is pinned to
the upstream Symbolica-3-compatible revision because the published 0.1.7 API predates it. `docs-notebooks` stages generated assets
under `docs/generated/notebooks`, which `docs-site` and `docs-watch` include automatically.
Rerun the export after editing a notebook; the watcher reloads the changed assets. The output
argument can instead point to an already built site. Generated wheels and notebook assets
remain local build outputs.

For the shared Spenso + Idenso showcase, also supply the matching Linnet browser wheel for
network figures:

// docs-example: syntax
```sh
just docs-notebooks /path/to/symbolica-wasm.whl spenso docs/generated/notebooks /path/to/linnet-wasm.whl
just docs-notebooks /path/to/symbolica-wasm.whl idenso docs/generated/notebooks /path/to/linnet-wasm.whl
```

Publication builds both wheels at the documented revision and adds all three products' assets
to the versioned site. The source notebooks use the checkout's model fixtures; the browser
export bundles those same inputs so it needs no external model download.

== Verify an installed host

Build and install the host's combined wheel, then run the repository smoke program in that
environment:

// docs-example: syntax
```sh
python crates/feynkit-py/tests/installed_import_smoke.py
```

It checks imports, representative classes, kinematics, and tensor reduction against the installed
Symbolica kernel. The repository also tests FeynKit/Spenso interoperability in both import orders,
but an in-process registration test does not replace this installed-wheel check.
]

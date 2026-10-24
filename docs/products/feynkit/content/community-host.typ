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

The 3.0.1 dependency versions alone do not select the Community display APIs
used by this checkout. Pin Symbolica, Numerica and Graphica together to the
compatible public Symbolica revision below, alongside the Typst atom payload.
These are the patches used by both this repository's root manifest and
`examples/notebooks/symbolica-host/Cargo.toml`. Put them in the consuming host's
root manifest so they apply to FeynKit and Spenso transitively, and retain the
resolved Cargo lockfile with the host release:

// docs-example: syntax
```toml
[patch.crates-io]
symbolica = { git = "https://github.com/symbolica-dev/symbolica", branch = "main" }
numerica = { git = "https://github.com/symbolica-dev/symbolica", branch = "main" }
graphica = { git = "https://github.com/symbolica-dev/symbolica", branch = "main" }
symbolica-typst-plugin = { git = "https://github.com/symbolica-dev/symbolica-typst-plugin", rev = "57b7455b37c30f4f87353c0073e49f5225e937ab" }
```

Only the top-level workspace's patches apply; dependency workspaces do not
forward their patch tables. The root, notebook host and separate `tydenso/`
Typst WASM workspace each own their patch selection. When building a separate
Typst plugin, align its resolved kernel and payload revisions with the host in
that build root too, then rebuild both artifacts.
Changing a Cargo patch does not update an already compiled Typst plugin or
native renderer. The notebook MathML path does not use the plugin.

In the host's existing core module, after `create_symbolica_module(m)?`, use its registration
macro:

// docs-example: syntax
```rust
register_module!(m, feynkit_py::FeynkitModule);
```

The host creates `symbolica.community.feynkit_native` and supplies `initialize_module()`.
Copy the package containing
#source-link("crates/feynkit-py/python/symbolica/community/hepkit/__init__.py", label: "the Python wrapper")
into the host's `python/symbolica/community/hepkit` tree. Its wrapper imports the native module and
calls the initializer. FeynKit itself must not declare another PyO3 extension entry point.
Registration must not construct symbolic expressions: callers need to import Symbolica before
calling `set_license_key`. Keep expression-valued constants as static methods. After rebuilding
the native host, run `python crates/feynkit-py/tests/installed_import.py` against that installation
to check that importing `E` and `hepkit` defers the license banner until the first operation.
Notebook figures use native Rust SVG geometry and direct MathML labels measured
in the notebook browser by AnyWidget. The native extension embeds no Rust Typst
compiler. Complete expression HTML also emits direct MathML; `formatted()`,
`paged()` and expression `to_html()` reject custom `notation_source`.
Bare `to_typst()` remains source-only;
`to_typst_project(..., notation_source=...)` exports `main.typ`, `render.typ`,
`notation.typ` and `tree.cbor` as a `dict[str, bytes]` without a compiler.
Expression `to_svg()` and graph `view.outlined()` use optional Python typst-py
(distribution `typst`).
The host supplies the verified font resources
listed in `crates/linnet-render-py/data.json` through
`SYMBOLICA_DATA_DIRECTORY` before rendering. This directory contains the extracted
`render/fonts/` assets; the renderer checks their compiled hashes.
Register `linnet_render_py::GraphModule` to install the shared
`symbolica.community.graph_native` module from `linnet-render-py`. Package
`crates/linnet-render-py/python/symbolica/community/graph` alongside HepKit so
`symbolica.community.graph` exposes its shared settings, views, snapshots and type information.

GammaLoop and the isolated notebook host in this checkout declare optional
`notebook-display` (`anywidget>=0.9,<1`), `graph-export`
(`playwright>=1.55,<2`) and `typst-export` (`typst>=0.15`) extras. These extras
are not currently released in the external Community assembly. Install the
equivalent dependencies when integrating another host; headless export also
requires installing Chromium:

// docs-example: syntax
```sh
python -m pip install "anywidget>=0.9,<1" "playwright>=1.55,<2" "typst>=0.15"
python -m playwright install chromium
```

`owner.render()` returns a `DiagramView` for generic and physics graphs,
networks, processes, amplitudes, generation collections, CFF and LTD. Displaying
it measures labels in the actual browser, then calls the existing native
`Scene::render`. After completion, `view.measured()` returns an immutable
snapshot whose `to_svg()`, `to_html()` and SVG/HTML `save()` retain the measured
boxes and embedded font. Scripts opt in with `view.measured(browser="chromium")`,
or obtain `view.measurement_request()` and return that same view's receipt to
`view.measured(measurements=receipt)`. A fresh owner's `to_html(browser=None)`
accepts no receipt and raises while measurement is pending; no browser is
launched implicitly and no estimated label bounds are substituted.
Default graph SVG contains XHTML `foreignObject` MathML for browser use.
PDF, PNG and Typst SVG import require `view.outlined()` so labels become vectors.

Native consumers prepare scenes with `PyDiagramView::from_scene` or
`from_scenes`; the former synchronous `PyDiagramRender::from_scene` entry point
is replaced by this measurement lifecycle. `PyDiagramRender` is the completed
snapshot. Views cloned in Rust share completion state, and a receipt belongs to
the exact view that requested it.

Bindings call the tensor owner's borrowed `display::atom_to_html`,
`display::structured_to_html` and `display::format_atom` APIs directly.
They do not import canonical Python modules to locate native functions.
The complete Spenso `mathml::expression` export returns `Result<String,
ExportError>`; the bounded semantic `mathml::render` callback retains its
optional notation contract. Native label import and measurement validation
return `LabelError`, while process scene construction returns `DiagramError`.
Convert those failures to host diagnostics at the presentation boundary.

When adding another component or class extension, follow the
#link("https://github.com/symbolica-dev/symbolica-community/blob/cad32c8f589c8a29531414f94b3fbe58020e0308/docs/writing-extensions.md")[Community extension guide].
Use the actual owner types and one native host. Run the same interoperability
assertions in separate standalone and Community processes, including a standalone
host without canonical package aliases and with canonical imports blocked.

Include the UFO loader release documented in the #link("guides/showcases/ufo/")[UFO import example]
when offering raw UFO import. Install it with `--no-deps` into the existing host environment:
the host supplies the single Symbolica extension. The checkout's `uv.lock` instead installs a
standalone Symbolica 3 kernel for GammaLoop development; do not sync that environment over a
community host containing its own native extension.

== Export the documented Python surface

Forward `feynkit-py/python_stubgen` from the host's stub-generation feature, alongside its other
community modules. Generic graph conversion, graph-valued topology callbacks, and Linnet selection
interoperability additionally require `feynkit-py/linnet-python` and registration
of `linnet_py::GraphModule` in place of the display-only module. Leave that
feature off in a display-only Community host. Python drawing selectors also
require the generic module, since their callbacks receive live graph views.
The documentation exporter updates both the package stub and the reference
input from the same native module:

// docs-example: syntax
```sh
cargo run --locked -p alphal00p-docs-python-exporter --features feynkit -- \
  feynkit-community docs/api/python/feynkit-community.pyi
cargo run --locked -p alphal00p-docs-python-exporter --features feynkit -- \
  feynkit-community docs/api/python/feynkit-community.pyi --check
```

Export the display-only rendering stub:

// docs-example: syntax
```sh
cargo run --locked -p linnet-render-py --features python_stubgen --bin linnet-render-stubgen
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
The browser wheel uses Zstandard level 22 inside the wheel's ZIP archive. Installing
it requires Python 3.14 with Zstandard support, as provided by the pinned Pyodide
runtime. Its ABI3 tag describes extension compatibility, not archive decompression
support in older Python installers.

The wheel command uses the checkout's pinned Emscripten Rust toolchain and provisions the
matching Pyodide build environment through cibuildwheel. Then export the executable
cells and start the existing documentation watcher:

// docs-example: syntax
```sh
just notebook-wheel
just notebook-ufo-wheel
just docs-notebooks /path/to/symbolica-wasm.whl feynkit docs/generated/notebooks target/notebook-ufo-wheel/ufo_model_loader-1.0.1-py3-none-any.whl
just docs-watch feynkit 8117
```

Pass the actual Symbolica wheel filename produced by the build. `notebook-ufo-wheel` downloads
the released `ufo-model-loader==1.0.1` pure Python wheel, matching the notebook dependency.
`docs-notebooks` stages generated assets
under `docs/generated/notebooks`, which `docs-site` and `docs-watch` include automatically.
Rerun the export after editing a notebook; the watcher reloads the changed assets. The output
argument can instead point to an already built site. Generated wheels and notebook assets
remain local build outputs.

The shared Spenso + Idenso showcases use the same host wheel for network figures:

// docs-example: syntax
```sh
just docs-notebooks /path/to/symbolica-wasm.whl spenso docs/generated/notebooks
just docs-notebooks /path/to/symbolica-wasm.whl idenso docs/generated/notebooks
```

Publication builds the host wheel at the documented revision and adds all three products' assets
to the versioned site. The FeynKit notebooks use built-in model constructors, so neither
native execution nor browser export needs separate model JSON files. The UFO-loading
tutorial deliberately retains its raw UFO fixture to demonstrate importing external models.

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
== Citations

Call `symbolica.get_citations()` after a calculation to collect Symbolica's citation
and references from the community packages used in that process. Each entry is a
`Citation`; `entry.to_bibtex()` exports its bibliography record. Importing a package
alone does not count as use. The collection is cumulative, deduplicates references,
and does not reset usage. Individual HEP entry points do not expose separate getters.
]

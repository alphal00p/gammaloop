# Five-product documentation

This directory is the source for the GammaLoop, Linnet, Spenso, Idenso, and
Vakint documentation sites. The five sites share one schema, renderer, theme,
and publishing pipeline, but each product owns its manual, navigation, search
index, PDF, Rust reference, Python reference, and changelog.

## Authoring model

- Put long-form product material in `products/<product>/content/*.typ`. Every
  chapter exports one named content value; register that symbol, its durable
  route, navigation group, title, and search summary in
  `products/registry.toml`. Each product must provide a tutorial and at least
  two manual chapters.
- Keep package `README.md` and `CHANGELOG.md` files concise. The builder may
  import them as explicitly typed Markdown; they are never evaluated as Typst.
- Link to the product that owns a concept instead of copying its explanation.
  Linnet owns graph algorithms, Spenso owns tensor execution, and Idenso owns
  representation and algebra identities.
- Generated API material under `api/` is deterministic input. Regenerate it
  through the exporters and commit the result; do not hand-edit snapshots.
- Classify every note under `architecture/` in `developers.toml`. The explicit
  section, summary, document type, and status keep implemented designs,
  proposals, dated investigations, diagrams, and archived baselines visibly
  distinct. The builder rejects an unclassified architecture note.

The builder passes authored Typst to the standard Typst 0.15 CLI. It does not
link Typst's internal Rust crates. Typst emits one semantic HTML document per
registered chapter plus the complete PDF manual. The Rust builder adds the
responsive book shell, navigation, heading anchors, local table of contents,
search, and previous/next links. Paged colors and components live in
`products/shared.typ`; the web theme and behavior live in `assets/site.css` and
`assets/site.js`.

## Rust metadata

`alphal00p-docs-schema` defines the versioned neutral catalog. A catalog uses
an explicitly constructed, ordered `DocScope`; there is no process-global
inventory. Component exporters decide which supported entries are registered,
while the Rustdoc sidecars remain exhaustive for the configured crate and
feature matrix. Rustdoc always disables Cargo defaults and enables exactly the
features listed in the registry; any manifest default must therefore also be
listed explicitly, which keeps the rendered provenance complete.

The proc-macro crate is exposed as `alphal00p_docs`:

```rust
/// Evaluates one supported operation. This first sentence remains ordinary
/// Rustdoc.
///
/// = Details
/// The remaining comment is raw Typst markup in the generated catalog.
#[alphal00p_docs::func(format = "typst")]
pub fn evaluate(input: Input) -> Output {
    // ...
}
```

The opt-in attributes are `func`, `ty`, `trait_item`, `macro_item`, and
`scope`. They capture stable source identifiers, signatures, parameters,
returns, fenced examples, fields, variants, and trait members into static
descriptors. For Typst-formatted comments, conventional Rustdoc keeps only the
first prose sentence. Unannotated existing Rust comments remain Rust Markdown.
Receiver methods infer their owning type. A receiverless associated function
uses `#[alphal00p_docs::func(owner = "Type")]`; the generated descriptor checks
that owner against `Self`, keeping method source identifiers collision-free.

The production component exporters use annotated `scope` descriptors in the
internal catalog adapter and source-backed proc-macro attributes for every
supported Rust item. Each adapter names a real workspace source file and stable
item identifier; macro expansion reads that declaration's comments, signature,
parameters, return type, fields, variants, or trait members and records an
`include_bytes!` dependency edge. Exporters then insert those descriptors in an
explicit runtime-scope order. This exercises the same proc-macro/runtime-scope
path for all five products without adding unpublished documentation dependencies
to the independently publishable product crates.

## Python metadata

All five Python references are built from isolated
`pyo3-stub-gen::StubInfo` exporters and checked-in deterministic `.pyi`
snapshots. Each exporter runs separately; combining the PyO3 inventories can
activate incompatible feature and link combinations. A one-component adapter
then writes neutral JSON, and the renderer merges only those serialized JSON
catalogs. GammaLoop's exporter also registers its real `_gammaloop` objects in
a fresh module and compares every export, class member, callable parameter, and
visible default with the generated snapshot.

Linnet's package and documentation stubs are byte-identical outputs of one
StubInfo renderer, and its import test checks their `__all__` against the real
extension surface. Spenso, Idenso, and
Vakint document their registered `symbolica.community.*` modules and declared
opaque return types. These community modules require the external Symbolica
assembly and are not standalone Python distributions.

The checked `.pyi` files are inputs and optional downloads, not the visible
reference. The site renders each neutral Python catalog as structured module,
class, function, parameter, return, member, example, feature, and source
documentation under `reference/python/<component>/`.

## Commands

```text
just docs-check
just docs-site gammaloop
just docs-site all
just docs-site all snapshot v0.3.4
just docs-watch spenso
just docs-serve spenso
just doc
```

`docs-check` verifies generated Rust and Python inputs, the product registry,
and Linnet's real imported Python surface against its checked-in stub. It also
generates an acceptance module from all five neutral catalogs and every fenced
block in the authored Typst manuals: Rust examples are compiled, catalog
examples marked side-effect-free are executed, shell examples are parsed
without running their commands, and TOML blocks are parsed as data. Every
Python example is passed to Python's compiler without importing the documented
module, so community-module examples do not require a standalone Symbolica
package in the ordinary documentation job. Source Rustdoc blocks explicitly
marked `ignore` remain non-standalone explanatory fragments; every supported
item also carries a separate compiling acceptance example. Examples that start
Symbolica's licensed runtime are compile-only here; the product's provisioned
integration jobs exercise those runtime paths.
`docs-site` emits a ready-to-publish Pages artifact under
`target/alphal00p-docs` by default. `docs-watch` uses a Rust-native filesystem
watcher and local HTTP server; successful generations trigger browser reloads,
while a failed rebuild leaves the last good site available and continues
watching. Rust-backed changes refresh only the applicable isolated API
exporters. A GammaLoop binding change also runs its runtime/signature parity
check before publication, so a stale checked-in stub leaves the last good site
served instead of publishing stale API pages. On startup, checked-in `.pyi` snapshots are trusted watch inputs,
while the source-backed GammaLoop and Vakint generated references are refreshed
without Python. `docs-serve` is retained as a compatible name for `docs-watch`.
Typst rendering, filesystem watching, and HTTP/SSE serving are all Rust-native
and do not use Python. A Python interpreter is needed only when validating or
regenerating an actual PyO3 binding surface. The existing `just doc` command
remains the workspace Rustdoc build.

Every product bundle has stable chapter routes:

```text
products/<product>/<channel>/
  .note
  index.html
  tutorial/index.html
  manual/<chapter>/index.html
  reference/index.html
  reference/python/<component>/index.html
  reference/rust/index.html
  reference/rust/supported/<component>/index.html
  reference/rust/<crate>/...  # exhaustive Rustdoc sidecar
  manual.pdf
```

GammaLoop additionally publishes its generated Clap command tree and settings schema at
`reference/cli/`; Vakint publishes its compiled topology and dependency inventory at
`reference/topologies/`. The neutral catalogs are the primary human-facing API reference.
Checked `.pyi` files and exhaustive Rustdoc remain downloadable implementation artifacts rather
than substitutes for tutorials or manuals.

The same tree is used for the moving `latest` channel and immutable tagged
snapshots. All authored internal links are interpreted relative to the product
root and rebased by the builder for nested chapter routes.
The plain-text `.note` file is the stable deployment probe for each product and
channel; it records the product, channel, tag, commit, and canonical route.

Contributor-facing architecture is published separately from the project
manuals:

```text
developers/
  .note
  index.html
  search-index.json
  architecture/<note>/index.html
```

Markdown notes use the shared responsive shell, navigation, search, and table
of contents. Standalone architecture diagrams retain their purpose-built
layout inside a clearly labelled developer page. Source links are pinned to
the documented commit. Latest publication replaces only the `developers/`
tree and product `latest/` trees; historical product snapshots remain intact.

Normal CI runs deterministic pure Vakint examples. FORM-backed examples use
`just docs-vakint-form-check`; pySecDec comparisons are explicit or scheduled
integration jobs.

## Publishing and snapshots

`.github/workflows/docs-pages.yml` publishes all five products under
`products/<product>/latest/`. A recognized component tag instead writes the
complete workspace suite to `products/<product>/snapshots/<tag>/`, records the
tag, commit, and every component version, and rejects a non-byte-identical
rewrite of an existing snapshot. Publishing is serialized and preserves the
`gh-pages` history.

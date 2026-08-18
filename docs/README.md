# Documentation contributor entry point

The GammaLoop, Linnet, Spenso, Idenso, and Vakint documentation is authored in
Typst and built by the workspace documentation crates. This compatibility
README is an index, not a second copy of the manual. Read
[`CONTRIBUTING.typ`](../CONTRIBUTING.typ) before editing documentation.

## Where to edit

- Product chapters: `products/<product>/content/*.typ`
- Product routes and generated-reference inputs: `products/registry.toml`
- Shared Typst components: `products/shared.typ`
- Developer and architecture records: `architecture/*.typ`, classified in
  `developers.toml`
- Shared web presentation: `assets/site.css` and `assets/site.js`
- Generated API inputs: `api/` (regenerate them; do not hand-edit snapshots)
- Documentation design and delivery contract:
  `architecture/documentation-improvement-plan.typ`

New authored documentation must be `.typ`. Files named exactly `README.md` and
`AGENTS.md` are the only Markdown compatibility exceptions; do not add a
parallel `README.typ` or duplicate a Typst page in Markdown.

Each product chapter exports one named value and has one entry in
`products/registry.toml` with its stable route, title, summary, group, source,
and symbol. Link to the product that owns a concept instead of copying its
explanation: Linnet owns graph algorithms, Spenso owns tensor execution, and
Idenso owns representation and algebra identities.

## Build and check

```text
just docs-check
just docs-site gammaloop
just docs-site all
just docs-site all snapshot v0.3.4
just docs-watch gammaloop
```

`docs-check` validates registries, generated CLI/settings and topology data,
Python inventories, checked examples, source-format policy, links, and assets.
`docs-site <product>` produces an honest product-scoped preview. `docs-site
all` builds the complete portal and all five products. Output defaults to
`target/alphal00p-docs`; `docs-watch` continuously rebuilds one product and
keeps the last successful preview available after an error.

Ordinary pull requests retain a complete preview artifact and cannot update the
mutable `latest` site. While pull request 96 is open, pushes to its in-repository
`docs` staging branch are the explicit development exception: they publish the
complete site to `latest` so it can be reviewed at the public Pages URL. The
workflow checks that the pull request is still open, so the exception stops
publishing automatically when it is closed or merged. After merge, `main` is
the sole source for `latest`. Recognized release tags publish immutable
snapshots.

## Generated reference

The human-facing Rust, Python, CLI/settings, and Vakint topology references are
rendered from versioned, source-backed catalogs. Checked `.pyi` files, JSON
catalogs, and exhaustive Rustdoc are provenance or tooling artifacts, not
substitutes for manuals. Improve source docstrings and catalog metadata rather
than maintaining handwritten signature copies.

Canonical product routes have this shape:

```text
products/<product>/<channel>/
  index.html
  tutorial/index.html
  guides/<topic>/index.html
  reference/interfaces/index.html
  reference/index.html
  reference/python/<component>/index.html
  reference/rust/supported/<component>/index.html
  reference/rust/<crate>/...
  version-history/index.html
  version-history/<component>/index.html
  manual.pdf
```

GammaLoop additionally publishes `reference/cli/`; Vakint publishes
`reference/topologies/`. Authored task-oriented material belongs under
`guides/`, generated API material under `reference/`, and release notes under
`version-history/`. Contributor architecture is published under `developers/`
in a complete build.

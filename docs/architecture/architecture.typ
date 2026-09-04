= Architecture notes
<architecture-notes>
This directory contains contributor-facing implementation notes,
proposals, diagrams, performance investigations, and build-system
records. The public index at `/developers/` is generated from
`docs/developers.toml`; that registry gives every published note a
constrained lifecycle, owner reference, review/evidence metadata,
applicable products and topics, and review triggers. Current records can
also pin source symbols and content digests: a changed verified scope
requires review instead of silently publishing an old architectural
claim against a new build revision.

The build revision, review revision, and historical evidence revision
are separate facts. `current` and `proposal` records age against their
review dates; `investigation` and `archived` records retain their capture
context without pretending to describe HEAD. Superseded records keep a
route and replacement edge but are removed from primary navigation and
search.

The collaboration portal is authored in `docs/portal/`. One Typst bundle
emits the landing, About, People, Talks, Publications, and Citations routes.
Those sources load `docs/portal.toml`, `docs/talks.toml`,
`docs/products/registry.toml`, and the generated INSPIRE publication cache
directly. Page bodies use ordinary Typst headings, paragraphs, emphasis, links,
and data-driven loops. Their layout calls name semantic portal components rather
than CSS classes. Only `docs/portal/components.typ` maps those components, the site
shell, and interactive browser controls to semantic HTML and CSS selectors. The
Rust documentation builder validates those inputs, invokes Typst, copies assets,
federates search, and checks links; it does not own the public page prose or HTML
composition.

Project tutorials and manuals remain under
`/products/<project>/latest/`. Architecture notes may explain internal
ownership and execution flow, but they do not replace the supported
interfaces documented in those manuals.

All maintained prose in this directory is Typst. Do not add a Markdown
or HTML twin; migrate a legacy source into one canonical `.typ` record
and update its registry and inbound links in the same change.

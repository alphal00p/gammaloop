= Architecture notes
<architecture-notes>
This directory contains contributor-facing implementation notes,
proposals, diagrams, performance investigations, and build-system
records. The public index at `/developers/` is generated from
`docs/developers.toml`; that registry gives every published note a
constrained lifecycle, owner reference, review/evidence metadata,
applicable products and topics, and review triggers. Current records can
also pin source symbols and content digests: a changed verified scope
requires review before a normal site build or publication. During local
`just docs-watch`, the same drift produces a warning and the preview
continues to rebuild, so contributors can update the code and its
explanation together.

The build revision, review revision, and historical evidence revision
are separate facts. `current` and `proposal` records age against their
review dates; `investigation` and `archived` records retain their capture
context without pretending to describe HEAD. `just docs-check` reports
maintenance warnings for missing owners, overdue or missing reviews, missing
review triggers, and missing historical evidence. Site and watch builds still
validate registry structure, scope boundaries, and links; they leave those
non-blocking maintenance reminders to the explicit audit. Only verified-content
drift becomes a warning in watch mode. Superseded records keep a
route and replacement edge but are removed from primary navigation and
search.

Build provenance uses `@` in a Jujutsu workspace and `HEAD` in an ordinary
Git checkout. The Jujutsu query uses `--ignore-working-copy`, so a watched
build reads the recorded working-copy commit without snapshotting edits.
Explicit commit and timestamp environment overrides take precedence. Source
archives supply `ALPHAL00P_DOCS_GIT_COMMIT` and
`ALPHAL00P_DOCS_GIT_TIMESTAMP` (or `SOURCE_DATE_EPOCH`).
Each build generation captures one commit and its timestamp for rendering
and link validation. Commits or Jujutsu snapshots made during that build
take effect on the next generation.

A verified scope normally hashes its whole source file. An optional
`end_anchor` narrows it to the exact bytes from the start of `anchor` up to,
but excluding, `end_anchor`. Both markers must occur exactly once, and the
end must follow the start. Missing, ambiguous, or reversed boundaries remain
errors, including in watch mode. Unrelated edits outside a selected range
do not invalidate its digest; changes within it still require owner review
or an explicit no-impact attestation for `just docs-check` and ordinary
builds. The registry lists the selected ranges, while broader review
triggers identify related code worth reviewing; a clean digest alone does
not prove that every architectural dependency is unchanged.

The collaboration portal is authored in `docs/portal/`. One Typst bundle
emits the landing, About, People, Talks, and Publications routes.
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

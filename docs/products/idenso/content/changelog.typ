#import "../../shared.typ": release-note, source-link

#let changelog = [
= Releases and change history

#release-note([
  `crates/idenso/CHANGELOG.md` is the canonical history, but its newest recorded release is
  0.2.5 while the current crate manifest is 0.3.0. The documentation build must show this gap
  and must not infer a 0.3.0 release note from commit history or implementation comments.
])

Idenso migrations should pay particular attention to representation initialization, dummy-index
handling, canonicalization, and Dirac, metric, and color conventions. A rewrite change can alter
the shape or ordering of a symbolic result without changing the mathematical intent, so callers
that persist or compare printed expressions should verify their chosen canonical form.

Immutable documentation snapshots pair the changelog with the exact Idenso crate version,
repository revision, feature set, and generated Rust/Python catalogs. Until a 0.3.0 entry is
curated, the rendered historical page must describe the existing changelog as partial coverage.

Read the canonical source at
#source-link("crates/idenso/CHANGELOG.md", label: "crates/idenso/CHANGELOG.md").
]

#import "../../shared.typ": release-note, source-link

#let changelog = [
= Releases and change history

#release-note([
  `crates/spenso/CHANGELOG.md` is the canonical core history, but its newest recorded release is
  0.5.6 while the current `spenso` manifest is 0.6.0. The documentation build must expose this
  coverage gap and must not invent a 0.6.0 entry from commits. `spenso-macros`,
  `spenso-hep-lib`, and `spynso3` have independent versions and changelogs.
])

When upgrading, inspect the history for tensor structure, contraction, network parsing,
Symbolica integration, and Python-stub changes. Also check each companion component used by the
application: a core Spenso release does not imply a release of the macro, HEP, or Python
adapter package.

An immutable documentation snapshot records all component versions from the registry alongside
the repository revision. Until the 0.6.0 history is curated, the 0.5.6 changelog must be labeled
as partial coverage rather than presented as the complete history of the selected build.

Canonical sources: #source-link("crates/spenso/CHANGELOG.md", label: "Spenso"),
#source-link("crates/spenso-macros/CHANGELOG.md", label: "macros"),
#source-link("crates/spenso-hep-lib/CHANGELOG.md", label: "HEP library"), and
#source-link("crates/spynso3/CHANGELOG.md", label: "Python adapter").
]

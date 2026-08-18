#import "../../shared.typ": callout, source-link

#let changelog = [
= Releases and change history

#callout("Version coverage", [
  This version of Spenso is 0.6.0, but the published changelog currently ends at 0.5.6. Release
  notes for 0.6.0 are not yet available, so the history below is incomplete for that version.
  `spenso-macros`, `spenso-hep-lib`, and `spynso3` have independent versions and changelogs.
])

When upgrading, inspect the history for tensor structure, contraction, network parsing,
Symbolica integration, and Python API or signature changes. Also check each companion component
used by the application: a core Spenso release does not imply a release of the macro, HEP, or
Python adapter package.

For reproducible upgrades, record the versions of Spenso and every companion component used by
the application. Do not assume that changes between 0.5.6 and 0.6.0 are fully represented in
the available notes.

Release notes: #source-link("crates/spenso/CHANGELOG.typ", label: "Spenso"),
#source-link("crates/spenso-macros/CHANGELOG.typ", label: "macros"),
#source-link("crates/spenso-hep-lib/CHANGELOG.typ", label: "HEP library"), and
#source-link("crates/spynso3/CHANGELOG.typ", label: "Python adapter").
]

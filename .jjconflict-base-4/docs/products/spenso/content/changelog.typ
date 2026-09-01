#import "../../shared.typ": release-note

#let changelog = [
= Version history

#release-note([
  This version of Spenso is 0.6.0, but the published changelog currently ends at 0.5.6. Release
  notes for 0.6.0 are not yet available, so the history below is incomplete for that version.
  `spenso-macros`, `spenso-hep-lib`, and `spynso3` have independent versions and changelogs.
])

== Available histories

- #link("version-history/spenso/")[Spenso core versions]
- #link("version-history/spenso-macros/")[Spenso macros versions]
- #link("version-history/spenso-hep-lib/")[Spenso HEP library versions]
- #link("version-history/spynso3/")[Spynso3 Python-adapter versions]

When upgrading, inspect the history for tensor structure, contraction, network parsing,
Symbolica integration, and Python API or signature changes. Also check each companion component
used by the application: a core Spenso release does not imply a release of the macro, HEP, or
Python adapter package.

== Upgrade checklist

Compare tensor structure and contraction behavior, network parsing, Symbolica integration,
feature flags, and Python signatures for every component used by the application.

== Reproducibility record

For reproducible upgrades, record the versions of Spenso and every companion component used by
the application. Do not assume that changes between 0.5.6 and 0.6.0 are fully represented in
the available notes.
]

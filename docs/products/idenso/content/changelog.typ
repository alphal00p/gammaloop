#import "../../shared.typ": callout, source-link

#let changelog = [
= Releases and change history

#callout("Version coverage", [
  This version of Idenso is 0.3.0, but the published changelog currently ends at 0.2.5. Release
  notes for 0.3.0 are not yet available, so the history below is incomplete for that version.
])

When upgrading Idenso, pay particular attention to representation initialization, dummy-index
handling, canonicalization, and Dirac, metric, and color conventions. A rewrite change can alter
the shape or ordering of a symbolic result without changing the mathematical intent, so callers
that persist or compare printed expressions should verify their chosen canonical form.

For reproducible upgrades, record the Idenso version and enabled features with the calculation.
Do not assume that changes between 0.2.5 and 0.3.0 are fully represented in the available notes.

Read the available release notes at
#source-link("crates/idenso/CHANGELOG.md", label: "crates/idenso/CHANGELOG.md").
]

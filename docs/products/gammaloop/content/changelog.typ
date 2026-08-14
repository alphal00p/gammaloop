#import "../../shared.typ": callout

#let changelog = [
= Releases and change history

GammaLoop versions are identified by repository tags together with the version of the Python
distribution and the `gammalooprs` Rust package. These components can evolve independently, so
record all of them when reporting a result or asking for support.

== Choosing the matching documentation

The `latest` channel follows current development. Documentation under `snapshots/<tag>` stays
with a tagged release. Use the snapshot matching your checkout whenever command options, state
files, or API signatures differ from `latest`.

== Before upgrading

Release notes currently live with repository tags and their linked pull requests rather than in
one consolidated changelog. Review them before upgrading, especially when you depend on
persisted states, CLI scripts, Python bindings, or numerical results. Keep a copy of valuable
state directories and re-run validation samples after changing versions.

#callout("Information to include in a report", [
  Include the repository tag or commit, the Python distribution version, the `gammalooprs`
  version, the run card, and any relevant settings overrides. This makes a calculation easier
  to reproduce and distinguishes a software change from a configuration change.
])
]

#import "../../shared.typ": release-note, source-link

#let changelog = [
= Version history

#release-note([
  GammaLoop does not yet have a consolidated changelog. Repository tags and their linked pull
  requests are the available release record. This page makes that gap explicit instead of
  presenting development documentation as a complete release history.
])

== Components and version sources

- Repository revision: #link("https://github.com/alphal00p/gammaloop/tags")[GammaLoop tags and releases].
- Python distribution: #source-link("pyproject.toml", label: "pyproject.toml").
- Core Rust library: #source-link("crates/gammalooprs/Cargo.toml", label: "gammalooprs package metadata").
- Rust/Python application API: #source-link("crates/gammaloop-api/Cargo.toml", label: "gammaloop-api package metadata").

These components can evolve independently. Record all of them when reporting a result or asking
for support.

== Choosing the matching documentation

The `latest` channel follows current development. Documentation under `snapshots/<tag>` stays
with a tagged release. Use the snapshot matching your checkout whenever command options, state
files, or API signatures differ from `latest`.

== Upgrade checklist

Release notes currently live with repository tags and their linked pull requests rather than in
one consolidated changelog. Review them before upgrading, especially when you depend on
persisted states, CLI scripts, Python bindings, or numerical results. Keep a copy of valuable
state directories and re-run validation samples after changing versions.

== Reproducibility record

Include the repository tag or commit, the Python distribution version, the `gammalooprs` version,
the run card, and any relevant settings overrides. This makes a calculation easier to reproduce
and distinguishes a software change from a configuration change.
]

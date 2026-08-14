#import "../../shared.typ": release-note, source-link

#let changelog = [
= Releases and change history

#release-note([
  GammaLoop does not currently have a canonical root `CHANGELOG.md`. Product releases are
  identified by the root Python version, the `gammalooprs` package version, and repository
  tags. This manual must not synthesize detailed historical entries from commit messages or
  reuse the changelogs of Linnet, Spenso, Idenso, or Vakint.
])

The documentation build may publish immutable snapshots for existing GammaLoop tags and a
moving `latest` channel. A snapshot records the repository revision and all component versions
from the registry. If release notes are curated later, they should describe user-visible CLI,
state, Rust API, Python API, and physics-result changes, including migration requirements.

Until a canonical changelog is added, consult the repository's tags and pull requests and
verify behavior against the current source. The short README remains the installation and
entry-point summary; architectural truth belongs in
#source-link("docs/architecture/architecture-current.md", label: "architecture-current.md").
]

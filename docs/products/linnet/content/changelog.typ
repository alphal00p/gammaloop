#import "../../shared.typ": release-note, source-link

#let changelog = [
= Releases and change history

#release-note([
  `crates/linnet/CHANGELOG.md` is the canonical Linnet history in this repository. Its latest
  recorded release matches the current `linnet` package version. The standalone `linnet-py`
  distribution has its own component version and currently has no separate changelog here.
])

The changelog records user-visible graph, parser, drawing, and algorithm changes. In particular,
users moving across releases should check index ordering, DOT serialization, subgraph behavior,
and feature-gated integration changes rather than assuming a serialized graph or incidental
ordering is stable.

Published documentation snapshots render the changelog from the same repository revision as
the API catalog. The moving `latest` channel may include an `Unreleased` section; an immutable
release channel must identify the selected crate and Python component versions separately.

Read the canonical source at
#source-link("crates/linnet/CHANGELOG.md", label: "crates/linnet/CHANGELOG.md").
]

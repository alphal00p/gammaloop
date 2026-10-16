#import "../../shared.typ": release-note

#let changelog = [
= Version history

#release-note([
  Linnet and Clinnet have canonical histories. The independently versioned Python `linnet`
  distribution does not yet have a separate changelog in this repository.
])

== Available histories

- #link("version-history/linnet/")[Linnet Rust library versions]
- #link("version-history/clinnet/")[Clinnet command-line renderer versions]
- Python `linnet`: no standalone history yet

The changelog records user-visible graph, parser, drawing, and algorithm changes. In particular,
users moving across releases should check index ordering, DOT serialization, subgraph behavior,
and feature-gated integration changes rather than assuming a serialized graph or incidental
ordering is stable.

== Choosing the matching documentation

The `latest` channel may describe unreleased behavior. Documentation under `snapshots/<tag>` is
fixed to a tagged repository revision. Choose a snapshot whose component information lists your
Rust crate version, and check the Python `linnet` version separately when using the Python bindings.

== Upgrade checklist

Review changes to index ordering, DOT serialization, subgraph behavior, and enabled features.
Re-validate serialized graphs and any code that depends on traversal or iteration order.

== Reproducibility record

Record the repository revision, Rust and Python `linnet` versions, enabled features, and the input
graph format. These distinguish a library change from a parser, serializer, or binding change.
]

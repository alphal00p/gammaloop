#import "../../shared.typ": callout, source-link

#let changelog = [
= Releases and change history

The #source-link("crates/linnet/CHANGELOG.md", label: "Linnet changelog") records changes to the
Rust crate. The standalone `linnet-py` distribution has its own version and does not currently
have a separate changelog in this repository.

The changelog records user-visible graph, parser, drawing, and algorithm changes. In particular,
users moving across releases should check index ordering, DOT serialization, subgraph behavior,
and feature-gated integration changes rather than assuming a serialized graph or incidental
ordering is stable.

== Choosing the matching documentation

The `latest` channel may describe unreleased behavior. Documentation under `snapshots/<tag>` is
fixed to a tagged repository revision. Choose a snapshot whose component information lists your
Rust crate version, and check the `linnet-py` version separately when using the Python bindings.

#callout("Before upgrading", [
  Review changes to index ordering, DOT serialization, subgraph behavior, and enabled features.
  Re-validate serialized graphs and any code that depends on traversal or iteration order.
])
]

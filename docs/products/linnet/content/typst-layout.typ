#import "/crates/linnest/typst/docs/manual.typ": layout-concepts, layout-reference, placement-concepts

#let typst-layout = [
= Layout functions

`layout` and `layout-sequence` are top-level functions exported by `lib.typ`, not members
of a public `layout` module. Their generated entries come from the shared `layout.typ`
implementation and use the public export names.

== Concepts

#placement-concepts
#layout-concepts

== Qualified aliases

The exported `layouts` module exposes the same implementations for callers that prefer
qualified access; these aliases link back to the canonical top-level entries below.

=== layouts.layout

Use #link("#layout")[`layout`] for the shared signature and parameter reference.

=== layouts.sequence

Use #link("#layout-sequence")[`layout-sequence`] for the shared multi-pass reference.

#layout-reference
]

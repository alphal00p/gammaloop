#import "/crates/linnest/typst/docs/manual.typ": drawing-concepts, drawing-reference

#let typst-drawing = [
= Drawing functions

`draw`, `edge-halves`, and `to-cetz-edge-halves` are top-level functions exported by
`lib.typ`. They turn an already laid-out graph into CeTZ content and expose the edge
geometry used by custom renderers.

== Concepts

#drawing-concepts

#drawing-reference
]

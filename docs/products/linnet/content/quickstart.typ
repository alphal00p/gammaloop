#import "../../shared.typ": boundary

#let quickstart = [
= Build your first Linnet graph

Linnet's half-edge graph model is available through three surfaces. Choose the language that owns
the graph in your workflow; the same concepts—paired internal half-edges, explicit boundaries,
typed identifiers, and validated topology—carry across them.

#boundary("Rust", [
  Use the published `linnet` crate when a Rust application owns graph construction, algorithms,
  or serialization. This is the most complete and thoroughly tested interface.

  #link("quickstart/rust/")[Use Linnet from Rust →]
])

#boundary("Python", [
  Use `linnet_py` for DOT-backed graph inspection and algorithms from Python. The binding is
  currently a source-built developer preview rather than a published wheel.

  #link("quickstart/python/")[Use Linnet from Python →]
])

#boundary("Typst", [
  Use Linnest when a Typst document owns graph construction, layout, and the final CeTZ drawing.
  The package is currently bundled with the repository rather than Typst Universe.

  #link("quickstart/typst/")[Use Linnet from Typst →]
])

Begin with Rust for application code, Python for an existing Python graph workflow, or Typst when
the desired result is an editable figure. The #link("reference/interfaces/")[interface guide]
explains where layout, rendering, and interchange belong.
]

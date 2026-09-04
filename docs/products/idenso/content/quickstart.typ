#import "../../shared.typ": boundary

#let quickstart = [
= Simplify your first tensor expression

Idenso supplies representation-aware identities through both Symbolica's Python community module
and a native Rust crate. Choose the interface that already owns your symbolic expressions; both
routes simplify a Minkowski metric contraction while preserving its one free index.

#boundary("Python", [
  Use the bundled community module for notebooks, exploratory algebra, and Python pipelines. This
  is the shortest supported installation path.

  #link("quickstart/python/")[Use Idenso from Python →]
])

#boundary("Rust", [
  Use the crate when a Rust application owns Symbolica expressions and should compose Idenso's
  metric, Dirac, color, or cooking passes directly.

  #link("quickstart/rust/")[Use Idenso from Rust →]
])

The #link("reference/interfaces/")[interface guide] maps both routes to the same identity families
and explains their feature and registration boundaries.
]

#import "../../shared.typ": boundary

#let quickstart = [
= Build your first Spenso tensor

Spenso is available as a native Rust library and as a community module bundled with Symbolica's
Python package. Both interfaces preserve representation, variance, and abstract-index structure
alongside tensor data.

#boundary("Python", [
  Start here for notebooks and exploratory calculations. The published Symbolica wheel includes
  the native Spenso community module and needs no GammaLoop checkout.

  #link("quickstart/python/")[Use Spenso from Python →]
])

#boundary("Rust", [
  Start here when a Rust application owns tensor storage, contractions, or network execution. The
  published `spenso` crate exposes the complete native API.

  #link("quickstart/rust/")[Use Spenso from Rust →]
])

Choose Python for the shortest installation and Rust for the most complete library surface. The
#link("reference/interfaces/")[interface guide] maps both routes to their generated reference.
]

#import "../../shared.typ": boundary
#import "/crates/kurvst/typst/docs/manual.typ": manual as kurvst-manual

#let kurvst = [
= Kurvst curve and path API

Kurvst is GammaLoop's Typst-facing WebAssembly layer for Bezier and Hobby paths,
arc-length trimming, sampled path patterns, parallel curves, and CeTZ conversion.
This page embeds the complete maintained Kurvst manual used by the drawing pipeline.

#boundary("Geometry rather than topology", [
  Kurvst transforms curve geometry. Linnet and Linnest own graph topology and layout;
  GammaLoop's drawing templates combine those results with Kurvst paths and styles.
])

#kurvst-manual
]

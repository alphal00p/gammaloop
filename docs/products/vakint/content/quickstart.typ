#import "../../shared.typ": boundary

#let quickstart = [
= Canonicalize your first Vakint integral

Vakint recognizes and canonicalizes single-scale vacuum integrals through both a Symbolica Python
community module and a native Rust crate. Both routes below deliberately disable evaluation
backends so the first success needs no FORM, MATAD, FMFT, or pySecDec installation.

#boundary("Python", [
  Use the community module for notebooks, exploratory matching, and Python algebra pipelines. The
  current Symbolica wheel bundles Vakint directly.

  #link("quickstart/python/")[Use Vakint from Python →]
])

#boundary("Rust", [
  Use the native crate when a Rust application owns integral expressions, matching policy, and
  later backend orchestration.

  #link("quickstart/rust/")[Use Vakint from Rust →]
])

Continue with the #link("guides/evaluation/")[matching and evaluation guide] only after the
canonical topology is understood and a backend policy is explicit.
]

#import "../../shared.typ": boundary

#let quickstart = [
= Generate your first FeynKit diagram

Both quickstarts load the repository's normalized scalar model and generate tree diagrams for
one scalar splitting into two. This checks the model-to-diagram API without a UFO installation,
FORM, or a GammaLoop state.

#boundary("Rust", [
  Use the workspace facade to load the model, configure generation, and inspect finalized
  diagrams. The example is compiled against the same revision as this manual.

  #link("quickstart/rust/")[Use FeynKit from Rust →]
])

#boundary("Python", [
  Use a Symbolica community host that includes FeynKit. Its expressions share the host's kernel,
  so model algebra, CFF results, and tensor reduction interoperate directly with Symbolica.

  #link("quickstart/python/")[Use FeynKit from Python →]
])

Continue with the #link("tutorial/")[generation tutorial] for physical process choices and the
#link("guides/notebooks/")[notebook guide] for figures. The scalar fixture demonstrates API
behavior; it is not a prediction for a measured decay rate or cross section.
]

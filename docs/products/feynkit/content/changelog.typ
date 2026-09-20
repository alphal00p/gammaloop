#import "../../shared.typ": release-note, source-link

#let changelog = [
= Version history

#release-note([
  FeynKit's focused crates currently carry version `0.1.0`. This manual describes the source
  revision recorded by the documentation build, not a claim that every matching crate or
  Symbolica wheel has been published with these APIs.
])

Record the Git revision, model fingerprint, generation options, Symbolica version, and tensor
reduction settings with results. Python installations additionally depend on which community
modules the host wheel registers. Rebuild that host and regenerate its stubs when changing the
native API; importing an older installed wheel does not exercise a newly edited checkout.

The current toolkit provides canonical model and graph ownership, deterministic generation,
shared CFF arenas, native rank-20 vacuum tensor projection, and notebook diagram rendering.
Use the #link("reference/interfaces/")[component references] to inspect the exact shipped surface.

GammaLoop states written before the unified FeynKit model/graph format require regeneration;
that application-level state version is separate from FeynKit crate versions. See
#source-link("docs/architecture/architecture-current.typ", label: "the persistence contract")
before replacing an existing GammaLoop binary or state.
]

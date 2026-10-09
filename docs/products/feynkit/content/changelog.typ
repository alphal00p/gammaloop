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

Python symbolic constants are now static methods: use `Symbols.dimension()`,
`Symbols.edge_momentum()`, `Kinematics.s()` (likewise `t()` and `u()`), and
`Kinematics.external_momentum()` / `loop_momentum()`. All expression-valued `Symbols` attributes
follow this convention, as does the tensor module's `Nc()`. Deferring their construction lets
callers set a license key after import.

The default loop momentum basis of a cross-section diagram no longer lets a cut edge carry part
of the beam momenta, such as `K + P1`: GammaLoop rescales loop momenta about the origin of this
basis, which must lie inside every cut energy surface. Regenerated cross sections can therefore
route differently from earlier exports, while saved states and DOT files keep their stored
`lmb_id` and explicit loop-momentum-basis overrides still take precedence. Diagram ids, names,
numerators and symmetry factors are unchanged.

Cross sections with identical beams, such as $g g$, $gamma gamma$ or $e^- e^-$, no longer come
out divided by the beam permutations. The CLI and Python defaults symmetrized cross-section final
states, which kept one forward graph per beam exchange without counting the graphs it stood
for. `symmetrize_final` now applies to amplitudes only: cross-section generation rejects it.
GammaLoop state version 9 refuses every state saved before this change, so regenerate them.
`symmetrize_initial` keeps one forward graph per exchange of two identical beams and multiplies
its overall factor by the number of graphs it stands for, so without left/right symmetrization
the factors add up to the unsymmetrized total; `symmetry_factor` is the automorphism count of
the kept graph. Python's `generate_cross_section` no longer takes `symmetrize_final` and, like
the CLI, applies no bridge filter by default. Left/right symmetrization merges a forward graph
with an identical mirror representative into one diagram whose factor counts both, so every
diagram id appears once.

The current toolkit provides canonical model and graph ownership, deterministic generation,
shared CFF arenas, native rank-20 vacuum tensor projection, and notebook diagram rendering.
Use the #link("reference/interfaces/")[component references] to inspect the exact shipped surface.

GammaLoop states written before the unified FeynKit model/graph format require regeneration;
that application-level state version is separate from FeynKit crate versions. See
#source-link("docs/architecture/architecture-current.typ", label: "the persistence contract")
before replacing an existing GammaLoop binary or state.
]

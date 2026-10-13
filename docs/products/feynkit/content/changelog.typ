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

Live generation pruning now uses `prune_skeletons`, `prune_dressings` and `select`
callbacks with typed views. Pruning remains hereditary; selection runs only on complete
assignments. Python fields are `Particle` objects, and views expire when their callback
returns. Explicit graph snapshots remain usable. Batched callbacks run on the calling
thread and stop immediately on cancellation or an exception.

Physics generation now separates an immutable `Request` from execution and output stages.
`Model.process(...)` defines external states; `.amplitude(loops=Loops.total(...))` and
`.cross_section(loops=Loops.forward(...))` require explicit finite loop ranges.
`Restrict`, `Orders`, `TopologyRules`, `Cuts`, and `Symmetrize` own discovery criteria;
`Run` owns threads, cancellation, and progress. Shared defaults are minimal and unrestricted,
without implicit GammaLoop CLI presets.

`Request.explain()` reports planning criteria and infeasible partitions. `count()` returns
exact integer counts and rational inverse-automorphism and signed physical-weight totals;
`generate()` returns compact prepared diagrams. `DiagramSet.build(NumeratorOptions(...))`
constructs numerators, with `Naming.none()` as the default. Explicit `Naming.named()` uses
the shared `GL` prefix, or a supplied prefix. Name selection and named momentum bases require
that naming demand; stable-ID selection and keyed bases compute IDs when needed.

`BuiltDiagramSet.select(Selection(...), run=run)`, `.group(Grouping(...))`, and `.into_amplitude(...)`
are explicit later stages. `Selection(select=[])` selects nothing; an omitted selection
retains all diagrams. Strings of exactly 32 hexadecimal digits denote IDs, so a short
hexadecimal name such as `D0` remains a name. Grouping preserves source indices and physical
factors, including when no names were requested.

Names and source indices are assigned during `build` and retained by later selection and
grouping. Grouped names can therefore differ from the former eager API, while diagram IDs
and physical contributions stay the same; grouping bookkeeping retains the build-stage
source indices.

Amplitude squared-order constraints require one common partner satisfying every bound.
The resolved relation survives selection, zero removal, grouping, and `into_amplitude()`;
squaring includes only allowed class pairs. Constructing a new `Amplitude` from a raw diagram
list has no request policy, so use `into_amplitude()` to retain that relation.

A diagram name is optional (`None` for unnamed diagrams); JSON and annotated DOT distinguish
it from an empty name. Rust's `FeynmanDiagram::id()` now returns a `Result`, computing and
caching the ID when requested. Generated aggregate numerators are constructed and cached on
first use from the stored local fragments; explicit or imported aggregate numerators are
retained. Serialization requests both the ID and aggregate.

Model fingerprints now preserve symbolic namespaces and function attributes while remaining
independent of symbol-registration order. Regenerate annotated diagram exports from their
original models when upgrading from the previous fingerprint encoding.

Cross-section cuts now use terminal-rooted partitions independent of which side is enumerated
first. Threshold finalization removes the initial-state attachment tree from both sides before
recomputing their crossing boundary, including loop-carrying crossings. Regenerated inventories
can differ from earlier exports.

Old GammaLoop DOT files using the previous `particle_id`/`vertex_rule_id` layout are
unsupported; re-export them using the current annotated format. Compact DOT identifies
particles by `particle` or `pdg`. Mass-only edges are unsupported, including a `mass`
attribute inherited from edge defaults: select a scalar model particle and configure its
model mass. A DOT `mass` attribute does not define a particle or override its model mass.

The default loop momentum basis of a cross-section diagram no longer lets a cut edge carry part
of the beam momenta, such as `K + P1`: GammaLoop rescales loop momenta about the origin of this
basis, which must lie inside every cut energy surface. Regenerated cross sections can therefore
route differently from earlier exports, while saved states and DOT files keep their stored
`lmb_id` and explicit loop-momentum-basis overrides still take precedence. Diagram ids, names,
numerators and symmetry factors are unchanged. GammaLoop refuses an integrand whose basis origin
lies on or outside a cut energy surface at the run's kinematics, whichever way the basis was
chosen, instead of losing or doubling sample points. The refusal also covers a cut that is
kinematically closed at those kinematics, which used to return meaningless values, and incoming
momenta given outside their centre-of-mass frame; its message names the cause.

Cross sections with massless beams, such as photons or gluons, no longer panic at their first
evaluation when the flux factor is on: the flux uses $m = 0$ for them. With the flux on, zero or
more than two incoming particles now return an error naming `general.disable_flux_factor`.

Cross sections with identical beams, such as $g g$, $gamma gamma$ or $e^- e^-$, no longer come
out divided by the beam permutations. The CLI and Python defaults symmetrized cross-section final
states, which kept one forward graph per beam exchange without counting the graphs it stood
for. `Symmetrize.FINAL` applies to amplitudes only: cross-section requests reject it.
GammaLoop state version 9 refuses every state saved before this change, so regenerate them.
`Symmetrize.INITIAL` keeps one forward graph per exchange of two identical beams and multiplies
its overall factor by the number of graphs it stands for, so without left/right symmetrization
the factors add up to the unsymmetrized total; `symmetry_factor` is the automorphism count of
the kept graph. Python's `process.cross_section(...)` rejects final-state symmetrization and, like
the CLI, applies no bridge filter by default. Left/right symmetrization merges a forward graph
with an identical mirror representative into one diagram whose factor counts both, so every
diagram id appears once.

GammaLoop's `import graphs` and Python `import_graphs` now read DOT through FeynKit, replacing
GammaLoop's own DOT parser. Plain compact DOT imports again: the model's Feynman rules fill in
every vertex factor, propagator numerator and external-state projector the input omits, and
`num` and `projector` attributes replace single factors. Such a factor must use the tensor indices
of the model's rules, `hedge(h, k)` for the k-th index on half-edge `h` counted from one; earlier
GammaLoop graphs wrote `hedge(h)` for the first one, and import now refuses a graph whose
numerator leaves an index other than external colour uncontracted, instead of failing, or
exhausting memory, at integrand generation. Compact cross-section DOT is accepted too: `is_cut` may mark one internal edge
instead of two legs, and without `final_state` every separating cut is physical until
`--process-spec` selects. `FeynmanDiagram.from_dot` resolves unqualified symbols in the
`gammalooprs` namespace, where the index, momentum and polarization heads live, instead of
creating unrelated `feynkit_graph::` symbols, and Rust's `DotDiagram` also returns the input's
attributes and which factors it wrote. Every edge needs `particle` or `pdg`; an edge's own `dir`
must agree with its particle; a vertex `int_id` must name a rule of the incident particles,
also next to `num`; a leg's `num` is ignored; and an edge `mass` replaces the particle mass in the
propagator denominator only. GammaLoop no longer reads `particle_id`, `vertex_rule_id`,
`initial_state_connection` or half-edge `ufo_order` payloads, and mass-only edges are gone.
`save dot` writes particle names and `int_id` instead, so its exports re-import; exports written
before this change name particles by `particle_id` and are refused with an error saying so, so
re-export them or replace `particle_id` by `particle`.

The current toolkit provides canonical model and graph ownership, deterministic generation,
shared CFF arenas, native rank-20 vacuum tensor projection, and notebook diagram rendering.
Use the #link("reference/interfaces/")[component references] to inspect the exact shipped surface.

Diagram ids hash a diagram's structural canonical form, with the names of the model's particles
and vertex rules: they no longer depend on model parameters, on the process that generated the
diagram or on its thread count, and they differ from the ids of earlier revisions. Select diagrams
again by the ids of the current revision.

Diagram validation, run by generation, numerator updates, cut and threshold installs and imports,
checks vertex slots and interaction signatures in one pass over the edges instead of one per
vertex, and each stored cut or threshold candidate in linear time. It performs the same checks
and reports the same first error as before.

GammaLoop states written before the unified FeynKit model/graph format require regeneration;
that application-level state version is separate from FeynKit crate versions. See
#source-link("docs/architecture/architecture-current.typ", label: "the persistence contract")
before replacing an existing GammaLoop binary or state.
]

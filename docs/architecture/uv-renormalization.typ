= UV renormalization architecture
<uv-renormalization-architecture>
#quote(block: true)[
#strong[Reviewed:] 2026-08-17 against `c9f4e32acd2c`

#strong[Lifecycle:] Current implementation architecture. Unsupported
prescription paths are recorded explicitly under
#link(<current-boundaries>)[Current Boundaries];.
]

== Scope
<scope>
This document records the implemented UV-renormalization invariants
shared by the legacy DAG-forest and hedge-poset orchestrators.
Historical revision comparisons, temporary log locations, and
external-reference convention tables belong in tests or investigation
records rather than in the current architecture.

== Orchestrators
<orchestrators>
`UVgenerationSettings::orchestrator` selects one of three execution
modes:

- `legacy_dag_forest` computes the established `Approximation` DAG and
  remains the reference implementation for connected graphs.
- `hedge_poset` is the default. It unfolds commuting UV operations into
  trace levels and stores their results in a compute store. It supports
  disconnected unions.
- `compare` computes both connected-graph backends, compares normalized
  expressions, and returns the legacy result. It is a validation mode,
  not a request to execute the returned integrand with the hedge-poset
  backend.

The backends share analytic operations and result types. Scheduling,
dependency lookup, caching, and disconnected-component composition
belong to the orchestrator because their traversal models differ.

The caller supplies the canonical model explicitly to every backend.
Forest construction, computation, and localization use that model for
particle, mass, and rule queries; the model does not own the traversal or
renormalization policy.

== Scheme and Computation Ownership
<scheme-and-computation-ownership>
The renormalization prescription belongs to the `Spinney`. A compute
node does not copy that policy; it stores computed values:

- the signed local four-dimensional counterterm;
- the integrated counterterm projections;
- cut-dependent local three-dimensional and final integrands.

This separation lets the same integrated result be consumed differently
without storing scheme-specific copies.

== Signs and Integrated Projections
<signs-and-integrated-projections>
Each UV operation supplies its own subtraction sign. In particular, the
local four-dimensional limit constructs the next counterterm as a signed
`-T(...)` operation. Nested signs therefore arise from recursive
operation composition, not from terminal depth or parity.

For a connected counterterm, integration stores a Laurent series from
which two physical projections are obtained:

```text
pole projection:               negative powers of integral(local)
finite-counterterm projection: -nonnegative powers of integral(local)
```

The second projection includes ε⁰ and the retained positive ε powers,
not only the strictly finite coefficient. Those positive powers are
required when factorized component series multiply into the finite part
of a disconnected counterterm.

No extra terminal parity is applied. Adding a sign derived from
operation count would duplicate signs already supplied by the local
operations.

The projections are used as follows:

- a non-root `PolePart` dependency recurses through the completed
  integrated pole only;
- a root `PolePart` result combines its local term with the integrated
  pole;
- an `MUV` dependency combines its local term with the signed finite
  counterterm;
- terminal renormalization output selects the pole or finite projection
  from the source Spinney\'s prescription;
- final three-dimensional integrand assembly always localizes the signed
  finite counterterm, regardless of the terminal prescription.

The last distinction is essential: changing terminal output policy must
not remove the finite addback required when assembling an integrable 3D
expression.

== Disconnected Composition
<disconnected-composition>
Disconnected composition depends on which object is being combined.

=== Four-dimensional local terms
<four-dimensional-local-terms>
Each connected component first obtains its complete recursive
counterterm, including the projection selected by that component\'s
Spinney. The union\'s local counterterm is the product of those full
component counterterms.

=== Integrated terms
<integrated-terms>
Integrated counterterms are scalars and factorize over connected
components. The stored disconnected value provides the product of
component pole projections and the product of component
finite-counterterm projections.

The Vakint return boundary restores Lorentz representations around source
abstract indices in momenta and metrics before further tensor operations.
Existing typed slots keep their dimension and identity; generated integer
indices remain scoped to the current UV node. A source label such as
`hedge(5)` must return as `mink(dim,hedge(5))`, rather than as an untyped
momentum argument mixed with typed tensor factors.

After tensor projection, the known regulated Minkowski slot dimension is restored
to the source representation before backend Laurent expansion. Representation
metadata stays constant during that scalar algebra, preventing derivatives of
metric/index declarations. Scalar evanescent factors still expand, including
inside function coefficients; unrelated representation dimensions are preserved.

The aggregate is sufficient when every component uses the same terminal
prescription. At terminal use, the hedge-poset backend nevertheless
returns to the component nodes, selects each component\'s own pole or
finite projection, and only then multiplies them. This also handles
mixed prescriptions without storing scheme-specific variants of the
aggregate.

=== Three-dimensional local terms
<three-dimensional-local-terms>
Cut CFF structures need not factorize over disconnected components.
Multiplying complete per-component local results would also multiply the
common root and can duplicate or incorrectly separate intertwined CFF
structure.

The hedge-poset backend therefore replays component-local operation
paths from their common root. At a union it:

+ obtains the replay states for every component;
+ forms compatible combinations of their integrated prefixes;
+ joins those prefixes in the trace/Foata representation;
+ applies each component\'s local operations to the shared incoming
  state.

The construction operates over an arbitrary number of components and is
not special-cased for a two-component spectacles graph.

== Marker Representation
<marker-representation>
UV markers are diagnostic symbolic structure. Approximation,
integration, series, and truncation are distinct canonical function
heads applied to the history of current/given subgraphs. Spenso printing
supplies the compact `K[...]`, angled-bracket, `Σ(...)`, and truncation
notation while structured DOT uses the same underlying atoms. Typst
presentation combines Symbolica\'s Typst arithmetic mode with those
Spenso shorthands; truncation is emitted as `op("Tr")(…)`, which is
valid Typst math, while compact terminal output remains `Tr(…)`.

DOT `full_num` fields contain Typst math fragments because the drawing
template evaluates them in math mode. Per-graph renormalization `.typ`
files wrap the same fragment in a display-math document, while the
accompanying `.txt` files remain the round-trippable Symbolica
representation.

Subgraph labels must be created through the subgraph symbol factory so
every consumer receives the same symbol metadata and print behavior.

== R\* preparation with UV counterterms
<rstar-preparation>

The experimental `uv::rstar` module implements the graph-first preparation
from the unfolded form of Chetyrkin's Eq. (4.49). Its input is a connected,
1PI homogeneous massless graph, initially IR safe, and one selected internal
line. In schematic notation the counterterm is

$ Z_G = -K I_("rem") T_("soft";c,M) [R'_U(G)
  + sum_(gamma in H_c(G)) C^U_(G,gamma) B_gamma] G. $

`RStarPlan` reuses the existing UV spinney catalogue for the ordinary
preparation. The separate IR kernel enumerates strict ancestor chains of each
proper hard component containing the selected line. Its sign is
`-(-1)^forest_length`; its UV calls act from smaller supports to larger
supports on the current quotient. Quotient connectivity and 1PI status are
tested after contraction. The full hard preparation includes its overall
counterterm, which supplies the mixed UV/IR contribution.

Every counterterm call retains its support, contracted insertion slots, open
boundary labels, degree bound, and projection scope. Ordinary child poles
belong inside the parent preparation's pole projection. In a quotient
c-operation, an insertion monomial is evaluated separately, with its Laurent
coefficient and hard loop variables outside that quotient's projection.

Lower-loop counterterms continue to use the ordinary massive R operation.
`RStarPlan.evaluate_uv_children` excludes the N-loop root before constructing
the old forest and retains the selected reference mass in the child
polynomials. `MassiveUvEvaluator` applies the same forest to a decorated
quotient, preserving the existing tensor insertion and contraction machinery.
The shared Vakint bridge unfolds traces and chains before legacy sum-interface
validation, so zero antisymmetric colour factors are removed before their
apparent tensor slots are compared. The R forest, Taylor expansion, and tensor
insertion algorithms are unchanged. The added forest accessor returns the
counterterm before any enclosing global projector.

At the Vakint handoff, GammaLoop identifies its Lorentz slot dimension `dim`
with `4 - 2 epsilon` before tensor reduction, including slots on opaque source
tensors. RustRed keeps these representation labels fixed while Laurent-expanding
scalar coefficients; GammaLoop restores its slot notation after integration.
Scalar evanescent factors still participate in the pole projection.

`TaylorExpansion` operates on an explicit expression whose boundary momenta
already carry a common scaling parameter. The outer substitution is
`mu_squared -> (1-lambda^2)*M_squared`; a remaining loop momentum in
`q + lambda*p` stays unscaled. The hard Taylor operation keeps its reference
mass fixed. Partial Taylor polynomials retain their scaling powers until their
product is truncated to the parent's total degree. A child may scale inherited
mass-squared parameters by `lambda^2` before introducing its own auxiliary
mass. Inserted polynomials participate in these substitutions even when the
selected line has been contracted.

`integrate_single_mass` first uses Vakint's existing vacuum tensor reducer,
then cuts the sole massive propagator. Vakint reduces the resulting massless
p-integral; `p_integral::MassiveClosure` performs the final radial integral
before expanding its epsilon-dependent momentum power. The same configured
normalization applies to every loop, including the closure. The GammaLoop
`forcer` evaluation-method setting selects the branch's p-integral backend.
Remaining massless vacuum integrals vanish, whereas already local terms are
retained. In particular, a swallowed selected line does not justify dropping
a mixed UV/IR polynomial or a quotient UV counterterm.

Each assembled remainder requests only its complete product's poles. An
insertion coefficient with epsilon valuation `v < 0` increases the scalar
integral's required depth by `-v`, preserving finite and higher terms that can
contribute to the final pole. Scalar hard-projector factors such as `1/d` are
included before this precision decision.

`RStarPlan.counterterm` assembles the preparation with general tensor
numerators and nonnegative superficial degree. Its root and hard Taylor
bounds count the total momentum degree, including externally forced vertex
momenta. For example, a ghost self-energy has loop UV degree one but needs
total Taylor degree two to retain its external antighost derivative times
the remaining linear expansion: `T_2[p F(p)] = p T_1[F(p)]`. The physical
vertex momentum is retained, and only the remaining factor needs a
first-order expansion. External projectors are excluded from this numerator.
UV classification still uses the original loop power counting. The existing vacuum tensor
projector first expresses the hard Taylor polynomial in complete local
covariants, leaving scalar hard-loop integrals unevaluated. For example,
`(k.q)^2` becomes `k^2 q^2/d`; the scalar coefficient `k^2/d` stays outside
the quotient UV map, while the complete covariant `q^2` enters it. The same
typed tensor algebra separates these factors for bare hard terms and local
counterterms. Between successive quotient maps, only scalar epsilon powers
are extracted; the complete tensor polynomial, including its metrics,
enters the next map.

This order corrects a tensor-specific ambiguity in the unfolded summary.
Keeping unprojected hard tensor components outside a quotient's pole
projection, while inserting an already tensor-reduced local counterterm,
uses different covariant bases for the two contributions. Since
`K(d/epsilon) != d*K(1/epsilon)`, that can change the final simple pole even
when their integrated hard preparation is UV finite. The hard vacuum
projection makes their bases consistent without computing their scalar
integrals or changing the original R tensor kernels.

All lower-loop UV maps use the ordinary massive routine, while the remaining
single-mass integrals use the p-integral path. Backend availability and loop
limits are those of Vakint. This is an explicit Rust API; the existing
renormalization command still selects its established R-operation implementation.

Exact regressions use the generated Standard Model Feynman rules. They compare
two selected-line choices for a two-loop gluon self-energy, and all eight
internal lines for a three-loop nested ghost self-energy, against ordinary
massive R counterterms. The
latter has a ghost line `A-B-C-D-E-F` and gluon chords `AF`, `BE`, `CD`,
testing ghost and gluon choices and multiple successive IR UV maps.
The one-loop ghost regression also checks the external-derivative Taylor
bookkeeping explicitly. No external projector is used in these comparisons.
With the fixture's `FMFTandMATAD` normalization and additional normalization
one, the ghost evaluations give
`GC_10^6 C_A^3 delta_ab (Q0.Q1) [-1/(48 epsilon^3) + 1/(32 epsilon^2) - 1/(48 epsilon)]`.
This retains the original external momentum orientations and omits global
graph prefactors, as does the counterterm API.

The four-loop nested ghost fixture extends the line to `A-B-C-D-E-F-G-H`
with gluon chords `AH`, `BG`, `CF`, `DE`. It compares all eleven internal
lines and the automatically selected line, and requires a
three-step IR chain. This is a selected-line consistency check; it has no
independent four-loop massive reference attached.
The counterterm is
`i GC_10^8 C_A^4 delta_ab (Q0.Q1) [-1/(384 epsilon^4) + 1/(128 epsilon^3) - 19/(1536 epsilon^2) + 5/(512 epsilon)]`
in the same fixture conventions.

The source identity is K. G. Chetyrkin,
#link("https://arxiv.org/pdf/1701.08627")[arXiv:1701.08627], Eq. (4.49)
in the PDF. Section 2.3 defines the smaller-first c-operation product with
earlier Laurent coefficients outside subsequent maps. Section 4.2 filters
each forest member by strict containment of the hard support; its nested
original supports label successive contracted pieces. The PDF identity is
derived for logarithmic graphs after differentiation. The compensated mass
Taylor expansion for power divergences follows the supplied _Rstar unfolded_
summary and is checked separately against the ordinary massive R operation.

=== Reconstructing the bare pole part

The local counterterm and the complete pole part of the original graph are
different when UV subdivergences are present. For an infrared-finite graph at
off-shell momentum, `RStarPlan.pole_part` uses
`K G = -Delta_U(G) - K sum_(S != empty) Delta_U(S) G/S`.
The sum runs over the same proper UV spinneys as the ordinary preparation.
Each counterterm retains its original tensor attachments. Its reference-mass
parameter is set to zero before insertion into the physical quotient.

Only the overall local counterterm uses the R\* rearrangement. The proper
quotients retain their original external momentum and massless denominators;
they receive no hard Taylor expansion. Their loop count is lower than the
original graph's, so the propagator backend can reconstruct five-loop poles
using at most four-loop p-integrals. Insertion poles are included before
Laurent truncation, retaining the required finite and positive-order terms
of those integrals. Tensor contraction also precedes the final pole projection.
An overall UV-convergent graph has zero local counterterm but can still have
nonlocal subdivergence poles, and is treated by the same reconstruction.

The API currently requires infrared-finite massless propagator graphs at
nonzero external virtuality and the lower-loop integrals supported by Vakint.
It does not recover physical IR poles of an IR-divergent input. Global graph
prefactors and external projectors remain outside both public methods.

The five-loop checks use entries 76, 159, 169, and 203 of
#link("https://arxiv.org/src/2104.08272v1/anc/MasterExpansions.m")[`MasterExpansions.m`],
the ancillary data of
#link("https://arxiv.org/abs/2104.08272")[_Glue-and-Cut at Five Loops_].
The fixtures preserve the published edge lists, have ten unit-power internal
propagators, and use the paper's G-scheme at `Q^2 = mu^2 = 1`. Each fixture
checks all ten internal massive lines plus automatic selection. GammaLoop's
`GScheme` setting selects Vakint's exact built-in normalization. All negative
epsilon powers are compared symbolically:

- Master 159: `36 zeta(3)^2 / (5 epsilon)`.
- Master 169: `441 zeta(7) / (40 epsilon)`.
- Master 76: `4 zeta(5)/epsilon^2 + (10 zeta(6) - 34 zeta(3)^2/5)/epsilon`.
- Master 203: `6 zeta(3)/(5 epsilon^3)
  + (9 zeta(4)/5 + 12 zeta(3)/5)/epsilon^2
  + (18 zeta(4)/5 + 96 zeta(3)/5 - 23 zeta(5)/5)/epsilon`.

The first two graphs are primitive. Master 76 has one proper UV bubble;
master 203 has two disjoint UV bubbles, exercising their product spinney
as well as the individual subtractions. Their proper UV children are one-loop
bubbles. Marking a line inside a bubble also exercises four-loop massive IR
quotient counterterms through Vakint's RustRed backend in symbolic master
mode. GammaLoop selects it with `evaluation_methods = ["rustred", "forcer"]`.
The engine explicitly loads the pinned four-loop external data; one through
three loops are vendored. One inside-bubble mark per nonprimitive master is
also checked using MATAD/FMFT as an independent reduction backend.
These exact reference comparisons use
`fmft.expand_masters = true` and `fmft.susbstitute_masters = false`: FMFT's
default substitution requests numerical constants. Symbolic FMFT and MATAD
output share exact zeta and harmonic-polylogarithm representations so these
constants cancel against each other and the p-integral result before the root
pole projection. The tests do not use a five-loop integral evaluation as
an oracle: the reference coefficients are the published ancillary values.
Lower-loop analytic checks retain the external scale logarithm, the first
open-index bubble moment, and the nonlocal pole of an overall-convergent
bubble-in-kite graph. Inverse powers of the external virtuality remain scalar
invariants during tensor restoration.

=== Choosing the massive line

`RStarPlan.with_automatic_edge(graph, model, settings)` chooses an internal
line before Taylor expansion or counterterm integration. It shares the original UV catalogue
between candidates and constructs each eligible hard component's IR forest
once. `RStarPlan.new` remains available for an explicitly chosen line.

The ranking first minimizes the largest edge-dependent loop reduction.
For each surviving remainder, it cuts the marked line, removes graph bridges,
and counts loops in the remaining connected components. This is a topology
estimate, not a guarantee that the backend will factor every numerator into
independent reductions. IR quotient UV calls
also contribute their loop counts. The undecorated proper UV children are
common across the choices and do not distinguish marks by loop count; their
mass-dependent polynomials and actual costs can still differ. This prefers a
cut that splits the propagator problem into smaller loop integrals, without
ignoring a costly massive IR quotient introduced by the same cut.

Ties minimize surviving graph summands: the ordinary UV terms plus, for each
hard component, its surviving preparation terms times its number of IR forest
terms. Further ties minimize quotient UV calls and then
degree-weighted hard preparations, then prefer a line outside every triangle.
The original edge index provides a deterministic final tie-break.
`estimated_cost` exposes these counts and the triangle flag. The cost
model and evaluator share the same scaleless test: an unintegrated loop
without the massive line vanishes, but a local mixed counterterm survives.

In particular, a line outside all proper UV components has no IR correction.
This can also avoid expensive massive quotient counterterms. Term counts are
structural estimates: numerator zeros, tensor moments, cache reuse and the
cut graph's propagator-reduction difficulty can change the actual runtime.
Preserving triangles retains potential local propagator reductions, but a
graph triangle alone does not certify a Forcer reduction rule: powers and
numerators also matter. The triangle preference breaks the otherwise equal
primitive-master scores. It picks edge 1 for master 159 and edge 7 for master
169, using the edge-list numbering in the fixtures. Counting how many
triangles an edge destroys does not reliably order the remaining choices,
so the tie-break only distinguishes cutting no triangle from cutting some.

IR ancestor enumeration also discards branches with a negative total local
degree before constructing their quotient operators. If the contracted
insertion has degree at most `d(gamma)`, the decorated ancestor's degree is
at most `d(A) - d(gamma) + degree(insertion) <= d(A)`. A negative bound
therefore makes that UV map zero. This uses the total local degree, including
forced external momenta and inherited mass powers, rather than membership in
the original loop-UV catalogue. The evaluator already rejected such
convergent decorated arguments; the earlier filter avoids their preparation.

== Current Boundaries
<current-boundaries>
- Local four-dimensional and integrated counterterm generation currently
  support only `MUV` and `PolePart`. The local three-dimensional kernel
  has an `IR` branch, but integrated `IR` generation is not implemented;
  `VaccuumLimit` and `OS` are also unsupported, while `Unsubtracted` is
  expected to be filtered out before these operations.
- Parametric integrand generation currently supports final 3D output
  only. `FourD` is used by integrated renormalization internally but is
  rejected by the parametric orchestrator.
- The legacy DAG executor constructs union-shaped nodes but does not
  execute multi-parent unions. Disconnected generation therefore
  requires the hedge-poset backend.
- `compare` cannot validate disconnected unions until the legacy backend
  has a corresponding execution path.
- External systems such as RQFT may use different propagator or vertex
  sign conventions. Those graph-specific conversion factors are recorded
  next to the reference expressions in renormalization tests; they are
  not part of the internal recursive sign convention.

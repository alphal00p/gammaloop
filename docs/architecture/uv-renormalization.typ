= UV renormalization architecture
<uv-renormalization-architecture>
#quote(block: true)[
#strong[Scope:] Shared soft/UV preparation with CFF and vacuum LTD projection.

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

== Physical loop measure
<physical-loop-measure>
The amplitude measure is $product_ell d^4 k_ell \/ \( 2 pi \)^4$. Signed CFF contours use $d q^0 \/ \( 2 pi i \)$, so each remaining CFF loop contributes $i \/ \( 2 pi \)^3$. Integrated UV subgraphs use the same Minkowski convention: `VakintSettings::additional_normalization` defaults to `"1"`; the historical `"-1"` inserted an extra minus sign per integrated loop. The standard Vakint normalization, forest subtraction signs, and causal prescriptions retain their separate meanings. A graph with $L$ loops consequently receives the same phase whether a loop is represented by CFF or by an integrated UV term.

== Dimensional numerator algebra
<dimensional-numerator-algebra>
The covariant `Local4dCts` representation retains symbolic Lorentz dimension `d`; its name does not mean that spin algebra has been evaluated at four. Analytic UV preparation reconstructs the subgraph's Feynman-rule numerator and changes its Minkowski slots to `d` before contracting metrics or evaluating ordinary closed fermion traces. The spinor trace normalization remains four. Expansion is permitted within the subgraph being analytically integrated, while the retained cograph and numerical-production numerator remain factorized.

GammaLoop owns the Dirac algebra. Vakint sees the Lorentz slots of opaque spin tensors and performs vacuum tensor projection in `d = 4 - 2 epsilon`. Projection can create new contractions between those slots. GammaLoop must restore and simplify the projected tensors in symbolic `d`, expose their scalar dimension dependence, and only then let Vakint select scalar backend orders and form Laurent coefficients. Merely collecting the resulting gamma chains does not evaluate their contractions. No hidden dimension-dependent trace may be deferred to final four-dimensional numerical assembly. The completed analytic result is expanded before its Laurent coefficients are stored and reinserted into the cograph, so equivalent scalar-master groupings produce the same counterterm. This normalization never expands the cograph.

Unresolved closed spin/Lorentz traces and internal contractions between open spin chains must report unsupported dimensional algebra. Retained free tensor slots remain a valid external basis. This avoids implicitly assigning a four-dimensional value to an evanescent operator; Sigma tensors also require a symbolic dimensional implementation before analytic UV use.

The required backend depth includes epsilon poles in numerator coefficients and loop-normalization factors. Positive epsilon orders needed by enclosing UV operations are retained using the forest's loop-order bound. A backend that cannot supply the required depth must report an error instead of silently returning a shallower expansion. The Laurent-series boundary also verifies integer powers and epsilon-independent coefficients. Supported scalar mathematical functions are expanded there; unresolved opaque epsilon dependence produces an error.

When integrated counterterms are enabled, the shared UV entry checks the raw edge and vertex numerators of the current UV subgraph before simplification. Gamma5 and chiral projectors are rejected with graph, subgraph and source context, including gamma5 pairs that could subsequently cancel. A dimensional gamma5 prescription is not implemented. Gamma5 belonging only to the cograph or an external-state projector is outside this check.

The generation setting `uv.project_integrated_uv_cts_onto_tensor_integrals` defaults to `true`. It selects Vakint's universal tensor kernels and isolated scalar coefficients. With `false`, each complete numerator follows the base branch's full FORM reduction and scalar-backend input path. Distinct denominator topologies remain separate; both modes use the same GammaLoop dimensional-algebra boundaries and conservative forest epsilon depth. Positive epsilon prefactors do not reduce that depth, and final truncation removes any excess orders.

== Orchestrators
<orchestrators>
`UVgenerationSettings::orchestrator` selects one of three execution
modes:

- `legacy_dag_forest` computes the established `Approximation` DAG and
  remains the reference implementation for connected graphs. It does not
  compute disconnected UV unions; those require `hedge_poset`.
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

Generation supplies all requested three-dimensional representations in
their configured order to one orchestrator call. Each backend unfolds its
forest once and computes the requested local four-dimensional Taylor terms and
integrated Vakint counterterms during the first representation pass.
Subsequent passes reuse those values and rebuild the cut-dependent 3D
projection and final integrand for the selected representation. Each
completed result is extracted before the next pass replaces the transient
cut store. `compare` performs this shared preparation independently for
each backend and compares the result for every requested representation.

One bounded projection cache remains local to the graph and backend across
all requested representations, so equal canonical Taylor sectors are reused.
Exact source and native-generation cache identities retain the full generation
options, including the representation and occurrence-local bounds. Each pass
reports only its own additional shared preparation time. Cached
payloads retain their certified numerator assignment; reuse never turns
independent Taylor sources into one capacity or coefficient sum. An LTD
source contributes its complete residue sum once, with its own numerator
maps, rather than once per possible production-orientation host.

== Scheme and Computation Ownership
<scheme-and-computation-ownership>
The renormalization prescription belongs to the `Spinney`. A compute
node does not copy that policy; it stores computed values:

- the signed local four-dimensional counterterm;
- the integrated counterterm projections;
- cut-dependent local three-dimensional and final integrands.

This separation lets the same integrated result be consumed differently
without storing scheme-specific copies.

Four-dimensional recursive state stores completed signed counterterm sectors
with their active component frames and frozen integrated frames. A containing
operator acts on each complete sector and extends the loop generators fixed by
its children. A spinney inherits the union of its components' loop generators;
its containing component supplies the remaining quotient directions. Surviving
child boundary momenta remain expressions in these coordinates and the current
component's external momenta. They are not additional required loop generators:
for a boundary momentum $q = r + p$, the containing soft operation scales
$p$ while retaining its loop coordinate $r$. Pure and nested terms must use a
consistent expansion basis; completed child loops cannot be redefined by an
external-momentum shift. No change of loop variables is applied to a completed
child to fit an enclosing basis.

The `has_soft_ancestry` flag records an
executed nonzero soft branch and propagates through later operators and
disconnected products; degree-zero H or a zero S branch introduces no flag.
The orchestrator validates scheme assignments before computation, while direct
`uv_limit` calls retain cheap invariant checks using this metadata. Ordinary
scheme controls are constructed explicitly in tests, without carrying a second
symbolic counterterm through production recursion.

== Nested soft-counterterm oracle
<nested-soft-toy>
The local-4D regression
`production_nested_toy_forest_terms_match_the_explicit_derivation` checks each
signed forest family before regrouping. The following scalar oracle makes the
physical-mass and completed-child dependence explicit. Use Minkowski products,
$V_ell = ell^2 - m^2$, $U_ell = ell^2 - M^2$, $r = k - p$,
$K = k dot.op p$, and $P = p dot.op q$, with child degree one and parent degree
two. Define

$ A = frac(k_0 + m, V_k V_r), quad C = frac(p_0, V_p V_(p + q)), quad X = A C. $

For $H = U + S - U S$, the child projections are

$ U_(gamma, 1) A = frac(k_0 + m, U_k^2) + frac(2 k_0 K, U_k^3), $
$ S_(gamma, 0) A = frac(k_0 + m, V_k^2), $
$ U_(gamma, 1) S_(gamma, 0) A = frac(k_0 + m, U_k^2). $

Thus the completed child is
$hat(A) = H_gamma A = (k_0 + m) / V_k^2 + 2 k_0 K / U_k^3$.
The soft part retains the physical mass $m$; $M$ enters through the UV part.
Write

$ A_0 = frac(k_0, U_k U_r), quad
  g_0 = frac(k_0, U_k^2) + frac(2 k_0 K, U_k^3), $
$ L_m = frac(1, V_p^2) - frac(2 P, V_p^3), quad
  Q_U = frac(4 P^2, U_p^4) - frac(q^2, U_p^3). $

Applying the parent to the bare and completed children gives the four families

$ F_nothing = frac(p_0 A, V_p V_(p + q)), quad
  F_gamma = -frac(p_0 hat(A), V_p V_(p + q)), $
$ F_Gamma = -p_0 A L_m - p_0 A_0 Q_U, quad
  F_(gamma; Gamma) = p_0 hat(A) L_m + p_0 g_0 Q_U. $

Their sum is

$ R = p_0 (A - hat(A)) (frac(1, V_p V_(p + q)) - L_m)
    + p_0 (g_0 - A_0) Q_U. $

The $g_0$ term is essential: the parent acts on the completed child, including
its UV contribution. Holding the child fixed would miss this nested term.
The test compares the four families independently and then checks their hard
and soft scaling; the regrouped sum alone cannot certify forest ownership.

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
- an `MUV` or `IR` dependency combines its local term with the signed finite
  counterterm;
- terminal renormalization output selects the pole or finite projection
  from the source Spinney\'s prescription;
- final three-dimensional integrand assembly always localizes the signed
  finite counterterm, regardless of the terminal prescription.

The last distinction is essential: changing terminal output policy must
not remove the finite addback required when assembling an integrable 3D
expression.

Vacuum masses carry hard weight within the component on which a Taylor operation acts. The four-dimensional route uses `mUV -> mUV/t` in its inverse-hard expansion; direct three-dimensional Taylor expansion uses `mUV -> lambda*mUV` with its forward hard scale. The auxiliary on-shell-energy deformation mass retains its separate role in the propagator rewrite.

A localized integrated contribution keeps its finite coefficient in the active expression. Its mass dependence scales when the current Taylor subgraph contains the coefficient's integrated owner, and stays fixed when the owner is disjoint. A later enclosing Taylor operation can therefore scale a coefficient that an earlier disjoint operation left fixed. Only the separate normalized localization kernel remains frozen outside all later Taylor operations.

The direct route retains this distinction with a transient one-argument `mUV(owner)` application of the existing vacuum-mass symbol. Each connected integrated coefficient is tagged before coefficients are multiplied. The owner uses the existing subgraph encoding: Taylor expansion scales the entire mass application for a contained owner, leaves a disjoint owner unchanged, and rejects partial overlaps. `DirectSector::combine` restores physical `mUV` on output copies; stored sectors retain their owners for further forest steps. This introduces no helper, type or symbol. Four-dimensional disconnected composition already applies the Taylor operations to their own components, so its runtime and `IntegratedCts` are unchanged.

Consumed loop measures are recorded separately. Writing `s = uvIntegratedLoopScale`, `IntegratedCts` projects a coefficient `C` whose integration consumed `r` loops as `C(mUV) * s^(4*r)`. The enclosing four-dimensional Taylor operation transforms `s -> s*t`, supplying `t^(4*r)` to compensate consumed loop measures without cancelling the coefficient's mass rescaling. Integration and physical projection set `s = 1`. The active denominator rewrite introduces the unscaled vacuum mass for the enclosing propagator basis.

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

The aggregate is sufficient when every component uses the same terminal
prescription. At terminal use, the hedge-poset backend nevertheless
returns to the component nodes, selects each component\'s own pole or
finite projection, and only then multiplies them. This also handles
mixed prescriptions without storing scheme-specific variants of the
aggregate.

=== Direct three-dimensional local terms
<direct-three-dimensional-local-terms>
Both direct local-3D representations perform the complete loop-energy integration first and apply every UV Taylor operator to the resulting complete/global CFF expression. They therefore share exactly the same local Taylor-transformed CFF bodies. If the generalized residue map is `{ k -> C_k }`, the orientation-parametric form is represented as `sum_k sigma(k) C_k`, where `sigma(k)` selects the complete residue-map key (loop map, edge map, numerator sampling map, and physical direction metadata). The Taylor operators treat `sigma(k)` as opaque. The sparse implementation therefore applies them branchwise without expanding the factorized numerator. `explicit_orientation_sum_only=true` merely replaces every `sigma(k)` by one and explicitly sums the same bodies. Neither direct form uses the exact-source reconstruction or minimax dispatch of completed local-4D Taylor terms.

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

Finite prefixes enter this replay as connected coefficients with their source
owners retained. Algebraic physical-mass factors carry that owner's mass tag;
scalar logarithms remain fixed coefficients of the soft jet. An enclosing
operation promotes the tags of contained children, including children whose
loops have already been integrated. A disjoint sibling with the same mass keeps
its own tag. The normalized frozen kernels remain outside the Taylor operation.
This preserves the engineering weight of finite mass vertices without rescaling
unrelated factors in a disconnected product.

The projected local-4D route is separate: it completes the Taylor expansion in four dimensions, reconstructs the exact source occurrence graph, and certifies the signed numerator mapping and denominator reconstruction before projecting that completed term to the requested representation. CFF selects a factorized EMR capacity assignment. Ordinary LTD retains the certified physical affine assignment without numerator-degree analysis; repeated LTD poles require bounds for exact finite sampling. Canonical denominator classes also gate the raw physical degree report, so diagnostic reporting does not reintroduce degree analysis for ordinary LTD. That reconstruction machinery is exclusive to projected local4D and is not a replacement for the direct replay above.

Final assembly retains the existing cut-keyed `Integrands` roots and shared
numerator definitions. Projected sectors complete their source-local residue sums
before the outer assembly; UV markers and tensor replacements use the existing
expression operations. CFF/from4D keeps its already-summed `SingleParametric`
construction. LTD uses the same generic evaluator builder for its complete sum,
with `SingleParametric` and `SummedFunctionMap` available when requested. The direct
CFF/from3D Taylor path keeps its existing keyed-family representation.
Within each direct Taylor operation, every maximal scalar subexpression that is
exactly homogeneous in the Taylor scale, such as the coefficient of a frozen soft
branch, becomes one shared closed factor times its power of that scale. No series
expands a factor body. The next enclosing operation inlines the factor and extracts
again in its own scale. The evaluator registers these factors as scalar functions,
not as numerator tensor components.

Requesting LTD requires local UV counterterms from expanded four-dimensional integrands and a complete orientation sum. Mixed CFF/LTD generation uses that same UV setting for both representations. Individual raised threshold residue-order coefficients may differ between CFF and LTD. Their complete contribution for each physical degenerate cut/residue must agree locally, with all pieces multiplying the same observable function grouped together. LTD threshold localization uses its own denominator factors and never reconstructs CFF to fix those individual coefficients.

== Finite-temperature zero-mode subtraction
<finite-temperature-zero-mode-subtraction>
At finite temperature a massless boson without a chemical potential keeps an
ungapped static Matsubara mode, which makes soft regions of the loop momenta
divergent. The local soft subtraction follows the effective-theory regions: it
keeps the static mode of the soft lines exact, Taylor-expands the gapped hard part
jointly in their momenta at fixed temperature, and adds nested zero bins where the
hard part itself contains a static line. Ownership and data flow are summarized in
#link("architecture-current.typ#3-2-three-dimensional-representations-and-numerator-energy-ownership")[the current architecture].
The host policy lives in `uv/thermal_zero_mode.rs`, the projection in
`uv/approx/direct_3d/zero_mode.rs`, the chart and node dispatch in `kernel.rs`, and
the frozen-branch tag rules in `thermal_soft.rs`. The examples below use the
three-loop ring graphs `tests/resources/graphs/thermal_*_ring.dot`: a quark bubble $A$
(loop $p$) and a second bubble $B$ (loop $k$: quark, gluon or ghost; a gluon tadpole
for the tadpole ring) joined by two static gluons carrying $q$; $h = A union B$,
and $G$ is the whole graph.

=== Soft sets, classes and orders
<zero-mode-soft-sets-classes-and-orders>
The soft set $Sigma(gamma)$ of a candidate host $gamma$ consists of the loop lines of
the cut subgraph outside $gamma$; bridges carry no loop momentum. It is eligible when
it is nonempty, every line is massless, has a zero mode (`UVE::zero_mode`: a boson,
ghosts included, without a chemical potential) and carries loop momentum but no
graph-external momentum, and it has $L_Sigma = L - L(gamma) >= 1$ loops, with $L$
the loop count of the cut subgraph. Its static infrared degree
(`UltravioletGraph::static_ir_degree`) is

$ omega(Sigma) = 3 L_Sigma + sum_(e in Sigma) "dod"(e) + sum_(v subset Sigma) "dod"(v), $

where the vertex sum runs over vertices whose every half-edge is a paired line of
$Sigma$, with the edge and vertex degrees of `local_dod`. The static mode diverges
iff $omega(Sigma) <= 0$.

Raised edge groups of $Sigma$ carry one Matsubara frequency at one on-shell energy
and form its frequency classes $c$ (`ZeroModeClassSpec`): $nu_c$ is the number of
edges and the smallest edge $r_c$ represents the class. Classes are taken by
decreasing $nu_c$, then $r_c$, while their signatures in a loop basis of the cut
subgraph with $gamma$ contracted stay independent, until there are $L_Sigma$. The
other lines of $Sigma$ are static conditions implied by the classes: they keep their
on-shell energies and receive no further temperature.

A host requests order tuples $a = (a_c)$ (`ZeroModeTuple`), each with the relative
Taylor order of the hard part

$ J(a) = -omega(Sigma) - sum_c (2 nu_c - a_c). $

- The static tuple $a_c = 2 nu_c$ reads the `ZeroModeSector::Complete` sector at
  $J = -omega(Sigma)$. A static mode with $J < 0$ converges and gets no host.
- Under vacuum subtraction a single class also has the odd vacuum-soft pole
  $a = 2 nu - 1$, read from the `ZeroModeSector::HardVacuum` sector, in which every
  distribution outside the class is at its vacuum value, at $J = -omega(Sigma) - 1$
  when that is nonnegative. The full thermal integrand is even in each class
  energy $x_c$, while the four-dimensional soft region of its vacuum part gives at
  most the pole $x_c^(1 - 2 nu_c)$, so tuples mixing static and odd orders vanish.
- `vacuum_subtraction = false` is supported: hosts then request only the static
  tuple.
- The projection is exact only for $a_c in {2 nu_c, 2 nu_c - 1}$. A host that needs
  a lower order (`UnsupportedOrders::BelowVacuumPole`) or the vacuum pole of several
  classes at once (`UnsupportedOrders::MultiClassVacuumPole`) is an error.

On the rings, the soft set of $h$ is the static gluon pair (one class, $nu = 2$),
with $omega = -1$: tuple $(4)$ at $J = 1$, plus $(3)$ at $J = 0$ under vacuum
subtraction. On the gluon and ghost rings the quark bubble $A$ also hosts
$Sigma(A)$, the static pair and the second bubble, with $omega = 0$ and the single
tuple $(4, 2)$ at $J = 0$.

=== Hosts and forest
<zero-mode-hosts-and-forest>
Candidates are connected `IR` spinneys of positive degree and disjoint unions whose
every component is one; `classified_spinneys` tags such unions `MUV`, so their
components decide. Hosts therefore need
`renormalization_prescription.massless_power_divergent = IR`. Per cut subgraph,
`assign_zero_mode_hosts`:

+ computes each candidate's canonical chart (`Graph::canonical_component_lmb`; its
  loop edges are stored as `coordinate_loop_edges`) and accepts candidates from the
  inside out. A line of the candidate with a static mode that carries crown momentum
  in this chart cannot be Taylor-expanded; it must lie in the soft set of an accepted
  inner host, its zero bin, or the candidate gets no region;
+ prunes an inner host contained in an outer host when its extra soft lines are safe
  in the outer chart, its static order does not exceed the outer one and it is no
  zero bin. Then $S^0_"out" S^0_"in" = S^0_"in"$ and its forest terms cancel; this
  removes the quark-bubble host of the tadpole ring;
+ requires the hosts to form a chain under inclusion, in which only a union host
  contains another, as one of its components, on a single uncut structure;
+ takes every cycle union, and the empty subgraph, as a hard part: an eligible soft
  set $Sigma'$ with $omega(Sigma') <= 0$ must be the soft set of a host, or contain
  the soft set of a host $gamma$ whose chart keeps the extra lines safe and with
  $omega(Sigma') + J_gamma + 1 > 0$, where $J_gamma$ is the host's static order;
+ installs the `ThermalZeroModeRegion` of a connected host on its spinney and that of
  a union host on an `IR` copy of the union (the plain union stays `MUV`), and
  demotes every other `IR` spinney to `MUV`.

The pass runs iff `subtract_uv` and
`UVgenerationSettings::freezes_thermal_soft_branches` hold: finite temperature,
`hedge_poset`, a 3D final integrand and `generate_integrated = false`. It does not
depend on vacuum subtraction. `Graph::restrict_full_observable_spinney` then keeps
the full-observable spinney under vacuum subtraction as a plain `MUV` U when a host
exists, and drops it otherwise; without vacuum subtraction the spinney keeps its
ordinary U. A host always needs it: its static cograph is logarithmically divergent
at the top Taylor order, and the whole-graph U is its local R\* subtraction. A host
without the retained spinney is an error.

A connected host and a union host $h$ apply

$ H_gamma = S^0_gamma + U_gamma - U_gamma S^0_gamma, quad H_h = S^0_h. $

The connected host's own U acts on its frozen terms in the
`Local3DDeformation::OrdinaryOnSoftBranch` role, with the branch scope that
enclosing operations co-scale, before the class propagators are multiplied in.

The union host is a second wood node (`WoodRole::UnionHost`) on the subgraph of its
plain factorized join, above it, and always one collective operation: every
transition into it carries its full scope. The U of each factor precedes it, and
the trace order applies $S^0_h$ to the outputs of the factors' U chains, which
supplies its UV subtraction. A factorized join replays its factors in a pinned
order: a factor holding a connected host whose soft set meets another factor
replays first, so the other factor's U acts on the zero mode ($U_B S^0_A$; the
reverse order gives $S^0_A U_B = 0$). Factors that must each precede another are
an error (`ZeroModeHostError::MutuallyCoveringJoin`). On the gluon and ghost rings
the host chains therefore sum to

$ (1 - S^0_h)(1 - U_B)(1 - H_A) F = (1 - S^0_h)(1 - U_B)(1 - U_A)(1 - S^0_A) F $

without any commutation assumption. The two-loop spinneys and the whole-graph U
complete the forest.

=== The zero-mode operation
<the-zero-mode-operation>
`apply_taylor` returns an all-zero input unchanged, as for an integrated prefix
covering one component of a union host. A nonzero input must expand the whole host
support in its stored canonical chart, with no soft line in the rescaled support;
violations are internal errors. For each tuple, `ZeroModeTerm::expand` then:

+ projects the tuple's sector. In the completed thermal CFF, the energies of all
  edges of class $c$ are set to one auxiliary energy $x_c$ in residue roots,
  numerator bodies and distribution arguments, and the class distributions become
  Laurent polynomials in $x_c$. The static mode is the pole $T \/ x_c^(2 nu_c)$
  times the static remainder; the other Matsubara modes are analytic at $x_c = 0$,
  since numerators and contact surfaces are polynomial in $x_c$. The projection is
  the Laurent coefficient
  $ P_(c, a) X = "OSE"(r_c)^(-a) [x_c^(-a)] X, $
  and $P_(c, 2 nu_c)$ selects the zero mode. Classes are projected in turn; the
  iterated coefficients are the joint Laurent coefficients. The temperature of the
  projected mode is written as the tag $T_("zm")(r_c)$. A summand of the form
  $f dot.op "OSE"(r_c)^(-a)$, with $f$ free of class energies, is read off
  unchanged, so the projection is idempotent; a frozen branch that would reach a `HardVacuum`
  projection unprojected is an internal error (`ZeroModeError::FrozenHardVacuum`);
+ Taylor-expands the coefficient, with the class propagators held apart, in the
  chart `Local3DDeformation::ZeroMode { endpoint }` with
  $"endpoint" = J - "dod"(gamma)$: the soft chart of the whole host support at fixed
  physical temperature and chemical potentials, through relative order $J$. It
  compensates every occupation of the support with one representative per
  statistics and chemical potential across all components of a union host, so the
  Landau pairs of the CFF variants cancel within each coefficient; the chart's series
  then cancels the momentum-dependent denominators of each coefficient symbolically,
  per numerator monomial (`NumeratorAtomExt::cancel_momentum_poles`), so the
  $1 \/ (p dot q)$ poles it creates never reach floating point. It co-scales the
  $T_"zm"$ of every class inside the support. A union host is one joint expansion of
  its components, a product truncated at total degree $J$, never a product of
  separate component expansions. On the rings the endpoints are $-3$ and $-4$ for
  $h$ (degree 4) and $-2$ for $A$ (degree 2);
+ re-keys every residue row to static soft energies
  (`DirectResidueBranches::rekey_static`), so numerators attached later see the soft
  lines at zero energy;
+ freezes the branch (`freeze_soft_branch`) and multiplies the class propagators
  back in as top-level factors, which an enclosing projection reads off without
  expanding the branch again.

=== Zero-mode tags and minimal depth
<zero-mode-tags-and-minimal-depth>
Two tags, defined next to `LOCAL_3D_SOFT_SCOPE` in `kernel.rs`, keep the static
structure visible until final assembly, which sets the temperature tag to $1 \/ beta$
and the order tag to one:

- `ZERO_MODE_TEMPERATURE` (`gammalooprs::zero_mode_temperature`), $T_"zm"$: the
  temperature of a static mode. Its argument is a class representative edge, as an
  integer, or the `S_<label>` owner of the completed zero-mode branch that took the
  class over.
- `ZERO_MODE_ORDER` (`gammalooprs::zero_mode_order`), $lambda_gamma$, with the
  `S_<label>` owner of $gamma$: the zero-mode chart writes it for the Taylor
  parameter of its finished series, so a term of relative order $j$ carries
  $lambda_gamma^(j - "dod"(gamma))$.

`normalize_local_tags` sets $lambda = 1$ and $T_"zm" = 1 \/ beta$ in the final
integrand; the overlap proof treats both tags as unfinished. Freezing a branch
re-owns the $T_"zm"$ of the classes inside its support, and those of strictly
nested zero-mode branches, to its own owner, and sets the strictly nested $lambda$
to one, so its own $lambda$ counts the whole vertex. A later operation multiplies a
tag by its Taylor parameter, like its loop momenta, according to the tag's role
(`SoftBranches::co_scales`; frozen means owned by a branch strictly inside the
operation's support):

#table(
  columns: (auto, 1fr, 1fr),
  table.header([Tag], [Zero-mode chart], [Ordinary charts]),
  [$T_("zm")(r)$], [co-scaled iff the class lies in the support], [fixed],
  [$T_("zm")(S_eta)$], [co-scaled iff frozen], [co-scaled iff frozen; an internal error on the own branch's U],
  [$lambda_eta$], [fixed], [co-scaled iff frozen; fixed on the own or a disjoint branch],
)

A frozen term $f_j q^j$ of relative order $j$, with its retained occupations and
carrier measure, is homogeneous of degree $"dod"(gamma) - j$ in every enclosing
chart under the co-scaling of frozen branches; the co-scaled
$lambda^(j - "dod"(gamma))$ compensates that degree, so only $q^j$ counts. An
enclosing U, whose series ends at order zero, thus expands unfrozen rows at their
four-dimensional depth and frozen rows at the UV degree of the static reduced
cograph, in one series call. On the rings: $U_G S^0_A$ has depth 0, $U_B S^0_A$
depth 1, the triangle spinneys give $U_D S^0_A = U_D U_A S^0_A = U_G U_D S^0_A = 0$,
and $U_G S^0_h$ expands a term of order $j$ at depth $j - 1$, so only the $q$-linear
term and the vacuum-soft term (depth 0) survive.

This minimal depth is a scheme choice, not a finiteness condition. It gives the
minimal add-back, the one of the dimensionally reduced effective theory, and the
depths of cancelling pairs coincide by construction. Expanding frozen rows at the
co-scaled four-dimensional depth would also leave a finite remainder, the
whole-graph U then covering the nonzero $U_D S^0_A$, but its add-back gains
$integral U_G U_D S^0_A$ and depth-two mass terms. Since a depth error stays UV
finite, UV profiles cannot detect it; the exact unit pins
`enclosing_u_acts_on_frozen_s0_at_the_static_cograph_degree`,
`minimal_depth_annihilates_subcritical_cographs` and
`vacuum_soft_term_gets_its_ring_u` fix the U images of the frozen branches. The first
and the last evaluate them with the static-theory mass below both equal to and apart
from $M$; the second checks symbolically that they vanish, whatever the masses.

=== Static-theory masses
<zero-mode-static-theory-masses>
The static soft lines of a residue row are the soft set $Sigma(gamma)$ of its
outermost completed zero-mode branch $gamma$: a branch strictly inside $gamma$ was
taken over by it, and the zero-mode operation of $gamma$ expanded that branch's
further soft lines as hard ones. A U whose support $X$ contains $Sigma(gamma)$ is, on
that row, the UV counterterm of the static theory: on the rings, $U_G$ on every row of
$S^0_A$ and $S^0_h$, the vacuum-soft rows $S^0_(h, 3)$ included. The masses it
introduces there, on the lines of $X$ that no inner U has made massive, are the
static-theory mass $m_"IR"$: the terminal symbol `mIR` and the expansion symbol
`mIRexp`, the roles of `mUV` and `mUVexp` (`OrdinaryEnergyMasses`), bound to the
runtime setting `general.m_ir`, which follows `general.m_uv` while unset. The choice
is made per forest term, not per U: the same $U_G$ keeps $M$ on $U_G F$ and
$U_G U_A F$. Every other mass is $M$:

- If $Sigma(gamma) subset.eq X$, $U_X F$ makes every line of $Sigma(gamma)$ massive,
  so it has no soft region of $gamma$, and $U_X S^0_gamma F$ cancels no infrared
  region of $U_X F$. With all its lines massive or gapped it is infrared finite on
  its own, and its UV cancellation against $S^0_gamma F$ holds for any mass: the
  rearrangement
  $1 \/ (l^2 + m^2) = sum_n (m_"IR"^2 - m^2)^n \/ (l^2 + m_"IR"^2)^(n + 1)$, truncated
  with the series, reproduces the large-momentum expansion through the depth of $U_X$
  whatever $m_"IR"$, at depth zero and one on the rings as at any higher depth, and so
  does an inner U with $M$ through its own depth.
- The own U of a connected host, $U_A S^0_A$ in $H_A$, keeps $M$, as the infrared
  requires: $Sigma(A)$ lies outside $A$, in the soft region of $Sigma(A)$ the term
  $U_A F$ tends to $[integral_p U_A Pi_A](q = 0)$ times the static mode, and
  $U_A S^0_A F$ cancels exactly that only with the same masses. A different mass
  leaves a constant times the static mode, a logarithmic infrared divergence on the
  gluon and ghost rings: with $m_"IR"$ there, $S(e_4)$ and its joint limits read
  about $0$ instead of $+1$.
- A U whose support holds only part of $Sigma(gamma)$, as $U_B S^0_A$ and
  $U_D S^0_A$, keeps $M$ as well, although no infrared cancellation fixes its mass.
  $U_X F$ has no soft region of $gamma$ either, its lines of $Sigma(gamma)$ in $X$
  being massive. The soft region that $U_X S^0_gamma F$ keeps, from its massless
  lines of $Sigma(gamma)$ outside $X$, is subtracted by the enclosing host acting on
  the same term ($S^0_h U_B S^0_A F$), and its UV cancellation against $S^0_gamma F$
  holds for any mass. Unlike the static theory's counterterm above, such a term is
  not infrared finite on its own, and $m_"IR"$ goes exactly to the terms that are.
  Giving it $m_"IR"$ instead also leaves every soft and joint limit and the UV
  profile of the ghost ring finite, with and without vacuum subtraction; only the own
  U taking $m_"IR"$ breaks them.
- Masses already in the row stay as they are: those inside the frozen branch, such as
  the $U_A$ of $[integral_p (1 - U_A) Pi_A^((0))]$, and those of an inner U on its soft
  lines, such as $U_B$ in $U_G U_B S^0_A$. An enclosing U rescales a terminal mass,
  `mUV` or `mIR`, like the loop momenta, and keeps its symbol.

$Sigma(gamma)$ consists of the loop lines of the uncut structure outside $gamma$
(`Graph::uncut_loop_lines`, on the only structure hosts live on), so a U containing a
frozen branch contains its soft set iff it contains every loop line; on the rings only
$U_G$ does. The ordinary U decides per residue row: `SoftBranches::row_owners` collects
the completed branches each row holds, visiting each definition its rows share once,
and `SoftBranches::expands_static_cograph` applies the rule. It expands the rows of a
static cograph in the chart `Local3DDeformation::OrdinaryOnStaticCograph`, the others
ordinarily, each part with its own copy of the numerator and scalar bodies its rows
call, so $U_G$ on $H_A F$ gives $m_"IR"$ to the rows of $S^0_A$ and $U_A S^0_A$ and $M$
to those of $U_A F$. A row holding a completed branch is a re-keyed static row of its
`S0`, apart from every row of the unfrozen integrand, so no row mixes the two; the
connected host, the one operation that adds both kinds of rows into one family, checks
that their residue keys stay apart. A row whose outermost branches have their soft
sets partly inside and partly outside $X$ is a generation error. The hosts form a chain
under inclusion, which leaves a single outermost branch per row, so no host assignment
produces one. The pins
`enclosing_u_acts_on_frozen_s0_at_the_static_cograph_degree`,
`vacuum_soft_term_gets_its_ring_u` and `vs_false_sector_mixes_full_and_minimal_depth`
evaluate every U image with $m_"IR" = M$ and $m_"IR" != M$;
`static_theory_mass_only_on_static_cographs` pins which images carry `mIR` and that
$U_G H_A F$ splits its rows, and `static_cograph_needs_the_whole_soft_set` pins the
row rule, the outermost branch deciding, and the error.

=== Integrated add-back
<zero-mode-integrated-add-back>
The integrated add-back of the zero-mode subtraction is not generated: with
`generate_integrated = true` the host policy does not run, and scheme validation
rejects every positive-degree finite-temperature `IR` component. The local remainder
$R$ needs no add-back to be integrable, but $integral R$ depends on the vacuum mass
$M$ (`mUV`) and the static-theory mass $m_"IR"$ (`mIR`) until the add-back is
supplied. On the rings, with $Pi_X$ the integrand
of bubble $X$ at static external $q$ after its Matsubara sum, $Pi_X^((i))$ its
order-$i$ Taylor term in $q$, $Pi_X^v$ its vacuum part, $N_0(k, q)$ the numerator
of the second bubble's vertices at zero energies, and $S^0_(h, 3)$ the vacuum-soft
tuple of $S^0_h$:

- Scaleless, with zero integral: $S^0_h$ and its $U_A$, $U_B$ variants
  ($integral q^(j - 4)$ and $integral |q|^(-3)$); $S^0_A$, whose UV and IR poles
  cancel (the two-loop logarithm of the static theory); $S^0_h S^0_A$,
  $U_A S^0_A$ and $U_B S^0_A$; and the odd whole-graph images $U_G S^0_h$ at
  $j = 1$, $U_G S^0_h S^0_A$ and $U_G S^0_h U_B S^0_A$.
- Nonzero and mass-dependent:
  - on the gluon and ghost rings,
    $ U_G (1 - U_A)(1 - U_B) S^0_A = & T^2 [integral_p (1 - U_A) Pi_A^((0))]
      integral_q integral_k \{ frac(N_0(k, q), (q^2 + m_"IR"^2)^2 (k^2 + m_"IR"^2)
      ((k - q)^2 + m_"IR"^2)) \
      & - frac(N_0(k, 0), (q^2 + m_"IR"^2)^2 (k^2 + M^2)^2) - (q"-linear") \}, $
    where the $U_A$ in the bracket and the $U_B$ behind the second and the
    $q$-linear terms keep $M$, with a $1 \/ epsilon$ pole and logarithms of $m_"IR"$
    and $M$, the matching pole of the static theory;
  - under vacuum subtraction, $U_G S^0_(h, 3) U_A U_B$, independent of $T$, with
    the $1 \/ epsilon$ pole of $integral (q^2 + m_"IR"^2)^(-3 \/ 2)$ and the vacuum
    vertices $U_A Pi_A^v (0)$, $U_B Pi_B^v (0)$ at $M$;
  - the ordinary vacuum add-backs of the U of $A$, $B$, the two-loop spinneys and
    $G$, at $M$.

No integration-level check of $d integral R \/ d M$ or $d integral R \/ d m_"IR"$
exists yet.

=== Limitations
<zero-mode-limitations>
- On the gluon and ghost rings the joint soft remainder ($q, k -> 0$ together)
  scales with power $+1$, from a term odd under $(q, k) -> (-q, -k)$; its
  parity-even part scales with $+2$. The integral converges, but with per-loop
  spherical maps the Monte Carlo variance of the odd term diverges
  logarithmically. Antithetic $(q, k) -> (-q, -k)$ sampling is not implemented.
- $omega(Sigma)$ ignores the suppression of a ghost vertex whose momentum factor is
  soft. It never misses a divergence, but on other graphs it can flag a convergent
  soft set containing ghost lines, which then gets an unnecessary but harmless
  counterterm or a spurious hostless error.
- A soft profile fits the successive differences of the remainder, so its $r^2$ gate
  fails on a ray where the remainder's linear approach coefficient vanishes by
  accident: the fit then mixes the linear and the quadratic order. The masses move
  these zeros. The static-theory terms are analytic where a hard quark line goes soft
  and only shift that coefficient, so $m_"IR"$ moves them as $M$ does: on the ghost
  ring with vacuum subtraction, the slow lane's $S(e_3)$ ray at $mu_B = 0$ reads
  $+3.95$ at $r^2 = 0.938$ with $m_"IR" = 0.3 M$, and without any static-theory mass
  (the code before `general.m_ir`) it reads $+3.84$ at $r^2 = 0.93$ for $M = 0.74$
  and $+4.08$ at $r^2 = 0.88$ for $M = 0.745$, while the remainder stays bounded. The
  slow lane's `*_with_static_theory_mass` cases therefore run at
  $m_"IR" = 0.5 M$, away from such a zero.
- A host needs the retained whole-graph spinney
  (`ZeroModeHostError::MissingEnclosingSpinney`), so graphs whose full-graph degree is
  negative are rejected: their static cograph has no local R\* subtraction.
- Near a soft static mode the remainder's rounding error comes from rounding the soft
  momenta against the hard ones and recurs at the rotated point, so the stability
  rotations cannot certify it. Evaluation therefore escalates up front at every point
  of a graph whose integrand subtracts a zero mode, that is, one where some counterterm
  hosts a zero-mode region; other finite-temperature graphs, bare ones included, keep
  the ordinary stability ladder. It starts in quadruple precision when
  the softest such line has an on-shell energy below
  `runtime.stability.zero_mode_energy_escalation_factor` times the temperature
  (default $10^(-2)$), or when the softest internal massless line of any species has
  one below `runtime.stability.zero_mode_hierarchy_quad_escalation_factor` times the
  sum $Sigma$ of the loop-momentum norms (default $5 dot 10^(-3)$). It starts at the
  last stability level when that line is below
  `runtime.stability.zero_mode_hierarchy_escalation_factor` times $Sigma$ (default
  $5 dot 10^(-6)$). The factors must be finite, and the quadruple-precision factor
  must be the larger one.
- These thresholds are measured on the four rings, with the static lines at an energy
  $E$ near $0.6 T$ and the default tolerance $10^(-5)$:
  - In double precision the remainder's error grows roughly like $(Sigma / E)^4$. With
    vacuum subtraction it first exceeds the tolerance at $E / Sigma$ between
    $1.4 dot 10^(-3)$ and $3.8 dot 10^(-3)$, depending on ring and ray, and stays below
    $2 dot 10^(-6)$ at $E / Sigma >= 5 dot 10^(-3)$.
  - A massless quark of energy $kappa < E$ multiplies that error by about $E / kappa$
    with vacuum subtraction and $(E / kappa)^2$ without it, without bound at a fixed
    $Sigma / E$: it is the bare integrand's own rounding loss, which grows like
    $kappa^(-2)$, amplified by the cancellation. Measured on the quark ring at
    $E / Sigma approx 6 dot 10^(-3)$, a quark at $kappa = 10^(-4) T$ already gives
    $10^(-4)$ with vacuum subtraction and $0.1$ without. That growth is slower than
    the growth in $Sigma / E$, so comparing the softest line of any species with $Sigma$
    keeps the error of every point left in double precision within that of a static
    line at the threshold.
  - Quadruple precision exceeds the tolerance from $E / Sigma approx 3 dot 10^(-6)$
    without vacuum subtraction and $3 dot 10^(-7)$ with it, and far earlier next to a
    soft quark.
  - Without vacuum subtraction, double precision still misses the tolerance by up to a
    factor of three just above the threshold, at $E / Sigma approx 6 dot 10^(-3)$, on a
    few rays with a hard second loop. The vacuum UV subtraction of the same graph,
    which has no static line, loses as much on those rays, so this is the rounding of
    the ordinary UV subtraction, which no a priori tier covers.
- Errors name the offending edge sets: a divergent soft set without a host
  (`Hostless`, naming `renormalization_prescription.massless_power_divergent = IR`),
  hosts that are not nested (`IncompatibleHosts`) or nested other than as a
  component of a union host (`UnsupportedNesting`), unsupported orders
  (`UnsupportedOrders`), and a host on a structure with cuts (`NonTrivialCut`).
  `CollectiveUnionScope`, `ClassRank` and the chart and tag checks of the kernel
  report internal errors.

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

== UV profiling from the CLI
<uv-profiling-from-the-cli>
`profile ultra-violet` probes the generated integrand along rays on which one or more loop momenta are scaled to infinity. It works for both amplitudes and Local-Unitarity cross sections. The ordinary, fast invocation is:

```text
profile ultra-violet -p <process> -i <integrand>
```

By default, `--selected-limits only-divergent` is active. For each graph, this mode uses the graph's expected bare UV degrees to retain only cycle unions with degree of divergence greater than or equal to zero. It profiles each distinct cycle union once in one generated physical loop-momentum basis, preferring the generation LMB when that basis represents the cycle. Thus the default avoids the old Cartesian scan over every loop subset in every LMB while still testing every expected divergent cycle. The target cycles are enumerated on the same physical domain as UV generation: dummy edges are excluded for amplitudes; LU first removes the incoming side and then enumerates separately in the complement of every compatible exact physical cut. If the stored physical LMBs cannot represent an expected divergent cycle exactly, profiling reports an error instead of silently claiming complete coverage.

The exhaustive behavior remains available explicitly:

```text
profile ultra-violet -p <process> -i <integrand> \
  --selected-limits all
```

`all` visits every nonempty loop-coordinate subset of every generated physical LMB. The same graph-theoretic cycle may consequently appear in several rows; that repetition is intentional because this mode audits all LMB representations. Cutkosky cuts do not create duplicate rows: in LU mode they are compatibility alternatives for a limit, not an additional profiling dimension.

One graph can be selected by its name, numeric ID, or `#`-prefixed numeric ID:

```text
profile ultra-violet -p <process> -i <integrand> --graph GL0
profile ultra-violet -p <process> -i <integrand> --graph '#2'
```

For a cross section, a particular physical Cutkosky cut can additionally be selected by its edge IDs. A graph selection is required:

```text
profile ultra-violet -p <process> -i <integrand> \
  --graph GL0 --cutkosky-cut 2,5
```

`--cut-edges` is a visible alias for `--cutkosky-cut`. The edge order is immaterial. GammaLoop rejects a missing, ambiguous, or amplitude-mode cut selection and reports the available physical cuts for a missing match.

=== Physical cuts and LU residue groups
<physical-cuts-and-lu-residue-groups>
Two identities must remain distinct when propagators are raised. Several physical Cutkosky cuts can correspond to the same normalized E-surface residue. Production LU evaluates that residue group once and emits an event tagged with the group's deterministic first cut. Therefore selecting any non-representative physical cut in that group correctly projects the weight of the representative event; evaluating another independent residue would double-count the group.

The selected physical cut nevertheless remains authoritative for deciding which UV limits exist. A cycle is compatible with a cut precisely when it does not contain an edge of that cut. Consequently:

- with `--cutkosky-cut` (or `--cut-edges`), compatibility is tested against exactly the requested physical edge set, even if another raised-edge alias belongs to the same residue group;
- without a cut selection, an LU limit is retained when it is compatible with at least one exact physical cut of the graph.

The second rule applies in both selected-limit modes. In exhaustive `all` mode it means that every nonempty LMB subset which is physically compatible with at least one cut is present, while the subset is still evaluated only once as part of the graph's summed LU integrand. In `only-divergent` mode the same compatibility filter is applied before retaining the unique expected-divergent cycle representatives.

=== Result and failure summaries
<result-and-failure-summaries>
The command first prints one detailed fitted-slope table per selected graph. The final `UV limit tests` line reports the passed and total counts. If any limit fails, it is followed by one colored rounded table containing only the failures, with columns for graph, selected cut edges, LMB, fixed and scaled loop edges, orientation, and failure reason. This final table is intended to make failures actionable without searching all successful rows. `--per-orientation` is available only for the orientation-parametric, localized direct-local3D representation; it includes orientation-local fits and failures in addition to the summed fit. Selector-free explicit-sum direct local3D and projected local4D are summed representations and reject this option rather than claiming per-production-orientation finiteness.

Each finite-precision ray is first fitted over its complete scale range. A missing or non-linear fit, or an invalid sample, triggers a complete Arb reevaluation of that ray. A sample that remains invalid prevents a fit and is reported as `missing_fit`; it cannot certify a zero tail. The profiler retains only the final weighted magnitude after evaluator-owned orientation summation and any selected-cut projection, and takes its logarithm before narrowing to the regression's `f64` inputs. The reported power is then fitted on the longest high-scale suffix with at least half of the samples (and at least five) whose $R^2$ is at least `0.99`; this permits a short pre-asymptotic transient without discarding isolated points. A fit which still has no such asymptotic suffix is reported as `unstable_fit`, even when its raw slope happens to lie below the UV threshold. At least five scale points are consequently required. The minimum and maximum scale exponents must also be finite and strictly increasing, so the fitted suffix is always the large-scale end of the ray.

Some rays vanish instead of following a power law. LU rays outside every compatible cut evaluate to exact zeros. The UV tails of vacuum-subtracted thermal integrands are suppressed by their occupations: their Arb samples fall faster than any power, then underflow to exact zeros. Cross sections, and amplitudes profiled with `--allow-vanishing-missing-fits`, accept a ray without an asymptotic power-law fit when it vanishes on a high-scale suffix at least as long as the fit window. A suffix with at most one nonzero sample must end in an exact zero; an isolated nonzero sample earlier in that suffix is allowed, but a nonzero endpoint is not evidence of vanishing. Otherwise, the nonzero samples must decay along their least concave majorant in the $(log lambda, log |f|)$ plane, with any exact zeros confined to the trailing samples. That majorant must fall faster than the threshold power from its first sample, and its last segment must be the final sampled step unless only exact zeros follow. For a nonzero tail, the initial envelope bounds at the start and at two later, evenly spaced sample indices must each become at least twice as negative. This conservative finite-window heuristic requires sustained acceleration; merely meeting the DOD threshold or having one steep final step is insufficient. It is evidence for thermal exponential suppression, not a proof of asymptotic behavior. The first bound describes decay over the sampled window, while the majorant bridges isolated dips. Such a dip occurs, for example, where Arb resolves an exact cancellation of a leading occupation at one scale and exposes the smaller remainder. A tail which flattens, rises or fluctuates above the threshold-power line is still reported as `unstable_fit`, and an asymptotic power-law fit always takes precedence.

For accepted tails, tables label the sparse terminal-zero evidence as `zero tail`. A finite envelope bound is labelled `decaying tail`, with `≤` and the measured bound in the slope column. Both leave the fitted DOD and $R^2$ blank rather than interpreting non-power-law evidence as a fitted power of negative infinity. With the exemption disabled, tables retain the failure and any available fit, noting the tail evidence without accepting it. The serialized `decay_slope` retains its finite bound or `"-inf"` marker for terminal-zero evidence.

== Current Boundaries
<current-boundaries>
- Local three-dimensional, local four-dimensional, and integrated counterterms
  support `MUV`, `PolePart`, and `IR` within the supported wood assignments.
  Positive-degree `IR` uses the completed soft prescription H; degree-zero H
  reduces to U. At finite temperature, positive-degree `IR` is the operation of a
  zero-mode host, described in
  #link(<finite-temperature-zero-mode-subtraction>)[Finite-temperature zero-mode subtraction],
  whose integrated add-back is not generated. Integrated `IR` uses the signed finite
  projection described above.
- Scheme validation rejects a positive-degree MUV parent containing a
  positive-degree IR child: its reduced-cograph subtraction lacks the containing
  soft branch. At finite temperature on the hedge-poset direct-3D route without
  integrated counterterms, the zero-mode host policy runs first: only zero-mode
  hosts keep `IR`, every other `IR` spinney becomes a plain U, and every enclosing
  operation treats a completed zero-mode branch as a frozen vertex, so a parent's U
  is the local subtraction of its static reduced cograph. On every other route a
  positive-degree finite-temperature `IR` component is rejected up front, including
  with integrated counterterms. It also rejects IR and PolePart components in one nested or
  disconnected UV wood until their combined integrated policy is implemented.
  Components that cannot occur together in a wood do not trigger this rejection.
- Local `OS` dispatch remains deferred, and requesting integrated OS is rejected
  before forest computation. `Unsubtracted` is expected to be filtered out
  before these operations. These checks are skipped when UV subtraction itself
  is disabled.
- Parametric integrand generation currently supports final 3D output
  only. `FourD` is used by integrated renormalization internally but is
  rejected by the parametric orchestrator.
- The parametric orchestrator receives the full generation settings so
  wood construction can see the medium. It rejects
  `local_uv_cts_from_expanded_4d_integrands` for medium modes and vacuum
  subtraction; their thermal-weight ownership is documented in
  #link("architecture-current.typ#3-2-three-dimensional-representations-and-numerator-energy-ownership")[the current architecture].
- LTD is vacuum-only and uses the projected local-4D route. Direct local-3D Taylor
  subtraction and orientation-local evaluation remain CFF-only.
- The legacy DAG executor constructs union-shaped nodes but does not
  execute multi-parent unions. Disconnected generation therefore
  requires the hedge-poset backend.
- `compare` cannot validate disconnected unions until the legacy backend
  has a corresponding execution path.
- External systems such as RQFT may use different propagator or vertex
  sign conventions. Those graph-specific conversion factors are recorded
  next to the reference expressions in renormalization tests; they are
  not part of the internal recursive sign convention.

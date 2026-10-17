= UV renormalization architecture
<uv-renormalization-architecture>
#quote(block: true)[
#strong[Reviewed:] 2026-10-01 against jj change `xuoqsyop`

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

== Symbolic values and backend boundaries
<symbolic-values-and-backend-boundaries>
Runtime edge and vertex numerators enter UV construction as
`SymbolicTensor<PartialStructure>`. `FourDSector::atom`, aggregate local-4D
values, `FourDTerm::numerator` and canonical projection numerators retain the
same carrier. Denominator classes, source incidence and momentum provenance
remain separate ownership data. Typed addition and multiplication preserve
logical ports and the rank of zero results.

Direct-3D `Integrands` keeps an Atom coefficient program beside shared
`FnMapEntry<SymbolicTensor<PartialStructure>>` definitions. The program's
numerator-family calls stand for those tensors; treating their opaque syntax
as evidence of a scalar interface would erase that distinction. Root
prefactors, markers and scalar energy substitutions remain separate from body
algebra. Final evaluator preparation combines the referenced interfaces before
lowering them through `SymbolicTensor::into_network` to finite components.
Only completed component functions become scalar `FnMapEntry<Atom>` entries.

Symbolica's series engine and the exact CFF energy mapper borrow Atom payloads.
The caller retains the source tensor and checks each returned coefficient or
mapped expression against its interface. Vakint is an explicit Atom backend
boundary: export decodes native index syntax, and returned projected tensors
receive checked admission before dimensional algebra resumes. Internal
materialization uses native `AbstractIndex` syntax; it does not decode to
GammaLoop's legacy index type and recook the result. `IntegratedCts` retains
typed pole and signed finite projections plus its loop-measure scale power.
`Series<AtomField>` is temporary calculation storage. Each completed projection
is checked against the source interface, including a vanishing projection.
Nested forests and the localizer consume those tensors directly.
These boundaries do not authorize distributing a graph numerator.

Symbolica seals `AtomCore`, so a downstream `SymbolicTensor` cannot implement
that trait and inherit its generic operation surface. Tensor-preserving
operations must be explicit in Idenso's API; broader trait interoperability
would require an upstream Symbolica change.

The rewriting APIs distinguish scalar algebra (`with_algebra_result`), checked
interface-preserving substitutions (`with_checked_expression`) and a genuinely
changed interface (`with_transformed_expression`). Matching ranks alone cannot
certify the identity of repeated anonymous ports. Custom tensor-valued heads
must declare their tensor semantics in addition to Symbolica linearity. Shared
definitions reuse their Arc only when both value and `same_reduction_state`
agree: a no-op contraction can establish a new completion certificate without
changing the expression. Typed stage handoffs keep those facts instead of
exporting and admitting the expression again. Timing evidence belongs in the
merge-readiness record; this ownership contract makes no speedup guarantee.

== Physical loop measure
<physical-loop-measure>
The amplitude measure is $product_ell d^4 k_ell \/ \( 2 pi \)^4$. Signed CFF contours use $d q^0 \/ \( 2 pi i \)$, so each remaining CFF loop contributes $i \/ \( 2 pi \)^3$. Integrated UV subgraphs use the same Minkowski convention: `VakintSettings::additional_normalization` defaults to `"1"`; the historical `"-1"` inserted an extra minus sign per integrated loop. The standard Vakint normalization, forest subtraction signs, and causal prescriptions retain their separate meanings. A graph with $L$ loops consequently receives the same phase whether a loop is represented by CFF or by an integrated UV term.

== Dimensional numerator algebra
<dimensional-numerator-algebra>
The covariant `Local4dCts` representation retains symbolic Lorentz dimension `d`; its name does not mean that spin algebra has been evaluated at four. Analytic UV preparation reconstructs the subgraph's Feynman-rule numerator and changes its Minkowski slots to `d` before contracting metrics or evaluating ordinary closed fermion traces. The spinor trace normalization remains four. Expansion is permitted within the subgraph being analytically integrated, while the retained cograph and numerical-production numerator remain factorized.

GammaLoop invokes Idenso's Dirac algebra before and after tensor projection. Vakint sees the Lorentz slots of opaque spin tensors and performs vacuum tensor projection in `d = 4 - 2 epsilon`. Projection can create new contractions between those slots. GammaLoop must restore and simplify the projected tensors in symbolic `d`, expose their scalar dimension dependence, and only then let Vakint select scalar backend orders and form Laurent coefficients. Merely collecting the resulting gamma chains does not evaluate their contractions. No hidden dimension-dependent trace may be deferred to final four-dimensional numerical assembly. The completed analytic result is expanded before its Laurent coefficients are stored and reinserted into the cograph, so equivalent scalar-master groupings produce the same counterterm. This normalization never expands the cograph.

Restoration reattaches the Minkowski representation to numeric backend labels and symbolic native labels. Only numeric labels acquire the topology-local index chart; native abstract labels retain their source identity. Already represented slots keep their complete interface.

Before reinsertion, completed internal momentum pairs are contracted and written
as scalar products. This prevents independent integrated factors from retaining
the same topology-local dummy label in their later product. Open boundary ports
remain available to the cograph, and its numerator is not distributed.

Unresolved closed spin/Lorentz traces and internal contractions between open spin chains must report unsupported dimensional algebra. Retained free tensor slots remain a valid external basis. This avoids implicitly assigning a four-dimensional value to an evanescent operator; Sigma tensors also require a symbolic dimensional implementation before analytic UV use.

The structural residual check materializes native indices with the existing
shorthand settings and strict slot-occurrence observer. Both encoded source
dummy names and native generated dummy names remain reserved. A genuine third
compatible occurrence is still an error; changing parser index types does not
relax multiplicity or unresolved-dimensional-algebra checks.

The required backend depth includes epsilon poles in numerator coefficients and loop-normalization factors. Positive epsilon orders needed by enclosing UV operations are retained using the forest's loop-order bound. A backend that cannot supply the required depth must report an error instead of silently returning a shallower expansion. The Laurent-series boundary also verifies integer powers and epsilon-independent coefficients. Supported scalar mathematical functions are expanded there; unresolved opaque epsilon dependence produces an error.

When integrated counterterms are enabled, the shared UV entry inspects the
source payloads of the retained edge and vertex tensors of the current UV
subgraph before simplification. Gamma5 and chiral projectors are rejected with graph, subgraph and source context, including gamma5 pairs that could subsequently cancel. A dimensional gamma5 prescription is not implemented. Gamma5 belonging only to the cograph or an external-state projector is outside this check.

The generation setting `uv.project_integrated_uv_cts_onto_tensor_integrals` defaults to `true`. It selects Vakint's universal tensor kernels and isolated scalar coefficients; `false` supplies the complete numerator to the selected tensor backend. This choice is independent of backend selection: the native FeynKit tensor reducer is the default, while AlphaLoop uses FORM. Scalar-master evaluation has its own backend policy. After tensor reduction, both projection modes use the same scalar coefficient kernel, retaining loop-independent factors outside the analytic backend and restoring them before final Laurent projection. Distinct denominator topologies remain separate; both projection modes use the same GammaLoop dimensional-algebra boundaries and conservative forest epsilon depth. Positive epsilon prefactors do not reduce that depth, and final truncation removes any excess orders.

Physical fixed momentum components are lowered to contractions with fixed basis
vectors at the Vakint boundary. Spatial basis contractions retain their Minkowski
signs; repeated `cind(0)` factors are independent temporal projections, not a dummy
index occurring more than twice. Open tensor ports and powered scalar spectators
remain owned by the shared tensor interface.

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

Comparison normalization sets the Lorentz dimension to four, reversibly cooks
index syntax for checked `SymbolicTensor` admission, reduces colour identities,
and contracts metrics and vector connections with chain and trace collection
disabled. It then converts scalar-product notation and restores the original
index encoding. Admission and normalization errors
propagate to the comparison caller. No Dirac or epsilon identity is enabled by
this normalization, and coefficient factorization is retained.

Exact finite-component certificates keep the tensor's external interface while
contracting its finite entries. Scalar normalization extracts numeric content
from integer powers before protecting their bodies. Functions, fractional
powers, and guarded inverse factors remain opaque. Positive scalar powers share
their complete base with linear occurrences; a summed base remains one opaque
factor. Restoring it may expose another linear cancellation, so the same
protected pass repeats only after a strict expression-size decrease.
Completed scalar components establish the absence of tensor interfaces once per
scope, avoiding repeated subtree scans. Tensor-bearing scopes keep the existing
factor-specific checks. The positive-power lookup skips summed-base
normalization when no such base was found; neither guard isolation nor the
collected expression changes.
The shared Spenso scalar collector owns literal numeric-content extraction;
GammaLoop uses that owner at its factor-preserving boundary.
Literal numeric content uses the real rational GCD of all real and imaginary
coefficient coordinates, including implicit units. The first coefficient in
coefficient-free body order chooses a phase in `{1, -1, i, -i}` so its remaining
real part is positive and its imaginary part is nonnegative. This preserves
opposite and quarter-turned power owners without invoking Gaussian-integer GCD
or changing any symbolic factor boundary.
Numerical weights are distributed once at each current product/sum boundary,
from the bottom up. Other multiplicative factors remain intact; rewritten
scalar subtrees are not recursively revisited. This avoids repeated copying
of deep signed sums in the scalar evaluator while preserving function,
inverse-power and guarded boundaries. Common numeric content is then
extracted by the same exact rule. Both enclosing normalization loops require
strictly smaller expressions before another pass.

Each finite-component pass makes one bounded collection in repeated literal
scalar factors. Complete scalar sums are admitted only after checking their
native scalar interface. A retained sum owns its proper additive subsets,
identified by their direct coefficient-free terms; selecting both would hide
their linear relation. Distinct summed candidates with the same
coefficient-free support also leave their constituents available to the
collector; complete positive-power bases remain opaque. Unrelated product
factors stay outside the exact coefficient rows. Recomposition can expose another cancellation, so the same
deterministic pass repeats only after a strict size decrease, at most eight
times. Unsupported collection or the pass limit leaves an unproven residual.
The certificate takes no numerical samples, forms no common denominators, and
never distributes a graph numerator.

Compact scalar collection also canonicalizes the precision metadata of an
exactly zero component in a complex floating coefficient to its nonzero
partner's precision. The nonzero component's value and precision bytes remain
unchanged. This permits exact reconstruction of expressions containing the
same numerical master coefficients without a tolerance or finite-value
rounding; it does not establish the analytic accuracy of those coefficients.
Source-frame certificates compare the independently normalized operands;
subtracting them first can itself recreate differing zero metadata.

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

For a connected counterterm, the native Laurent series is temporary calculation
storage. `IntegratedCts` retains two typed physical projections with the source
logical interface, including a vanishing projection:

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

The finite coefficient is stored as a shared numerator definition, using the existing residue-key owner to apply its component momentum map exactly once. Scalar CFF branch weights and the frozen kernel remain separate. A later cograph numerator is another shared factor; it does not absorb or remap the integrated tensor vertex. This permits open spinor coefficients to be added to local counterterms without asking scalar branch assembly to combine incompatible tensor interfaces. The localization regression checks both the scalar affine substitution and the exact product of separate finite and cograph matrices.

The direct route retains this distinction with a transient one-argument `mUV(owner)` application of the existing vacuum-mass symbol. Each connected integrated coefficient is tagged before coefficients are multiplied. The owner uses the existing subgraph encoding: Taylor expansion scales the entire mass application for a contained owner, leaves a disjoint owner unchanged, and rejects partial overlaps. `DirectSector::combine` restores physical `mUV` on output copies; stored sectors retain their owners for further forest steps. This introduces no helper, type or symbol. Four-dimensional disconnected composition already applies the Taylor operations to their own components; the retained-tensor change preserves that mass ownership and the projection signs.

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
Integrated coefficients have no remaining loop dependence from their
integrated components, but can retain open tensor indices. They factorize
over disjoint connected components. The stored disconnected value provides the product of
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

The projected local-4D route is separate: it completes the Taylor expansion in four dimensions, reconstructs the exact source occurrence graph, performs the factorized minimax EMR dispatch needed for derivative-created occurrences, and only then projects that completed term to CFF. That reconstruction machinery is exclusive to projected local4D and is not a replacement for the direct replay above.

== Marker Representation
<marker-representation>
UV markers are diagnostic symbolic structure. Approximation,
integration, series, and truncation are distinct canonical function
heads applied to the history of current/given subgraphs. Marker
calls remain actual arguments rather than hidden body constants, so an outer
approximation or cograph prefix updates each history before evaluator binding.
This does not distribute the numerator. Spenso printing
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
  reduces to U. Integrated `IR` uses the signed finite projection described above.
- Scheme validation rejects a positive-degree MUV parent containing a
  positive-degree IR child: its reduced-cograph subtraction lacks the containing
  soft branch. It also rejects IR and PolePart components in one nested or
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
  #link("architecture-current.typ#3-2-cff-production-and-numerator-energy-ownership")[the current architecture].
- The legacy DAG executor constructs union-shaped nodes but does not
  execute multi-parent unions. Disconnected generation therefore
  requires the hedge-poset backend.
- `compare` cannot validate disconnected unions until the legacy backend
  has a corresponding execution path.
- External systems such as RQFT may use different propagator or vertex
  sign conventions. Those graph-specific conversion factors are recorded
  next to the reference expressions in renormalization tests; they are
  not part of the internal recursive sign convention.

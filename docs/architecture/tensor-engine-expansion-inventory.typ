= Tensor engine expansion inventory
<tensor-engine-expansion-inventory>

This is the M0 inventory for R17 of
#link("tensor-engine-consolidation-plan.typ")[the tensor engine consolidation plan].
It records the fixed implementation base `bfb80d8c`, including the preserved
ladder benchmark and bounded sixteen-gamma trace work. It is an inventory of
work to do, not a claim that the internal expansion paths have been removed.
Line numbers below refer to that base; function names identify the owners
when subsequent changes move them.

== Current owner audit

The combined implementation was re-audited on 2026-09-28. The direct
production expansion calls in Idenso and Spenso now belong to explicit
materialization: `SymbolicTensor::expanded`, the alias materializer,
Spenso's component/network expansion operations, and filtered
`expand_with_map`. The scalar-sum emitter is called only by explicit
`SymbolicTensor::expanded`; contraction cannot dispatch to it. The flat
expanded collector and the trace-specific sparse emitter are test-only
oracles. Portable-payload `collect_symbol` methods record symbol metadata
and perform no algebraic collection.

Contraction reads the borrowed factor graph, visits one factor's alternatives
and merges equal remaining states while retaining alias definitions. Gamma
and colour passes retain scalar and foreign spectators. Collection shares
Spenso's term tape and keeps cofactors outside selected alternatives. The
legacy network Schoonschip engine and selective-expansion wrappers have been
removed. Physical adjunction parses the compact matrix structure without
expanding its body; the formal transposition mode alone requires a unique
external pairing across branches.

The historical site table below remains the removal inventory at its recorded
base. This source audit establishes ownership of the surviving operations;
the combined runtime, callback and component suites remain the functional
acceptance checks. It does not by itself certify plan completion.

== Boundary

Rule application may introduce a sum prescribed by an identity. Contraction
may visit alternatives of a factor and merge equal remaining states. Neither
operation may construct the distributed product of those alternatives. They
return alias definitions and a typed root. Only an explicit caller request
to materialize or expand may collect the complete polynomial.

The existing `SymbolicTensor::expanded` operation, component-tensor expansion,
and `Network`'s explicit `TensorAtomMaps` expansion methods are materialization
entry points. Their Symbolica calls remain valid when invoked explicitly;
engine passes must stop calling them. Component execution remains Spenso's
responsibility. Graph numerators retain factorization even in diagnostic or
test copies. Use factor-preserving symbolic certificates and independent exact
component checks there; only non-graph explicit materialization fixtures may
use expansion as their oracle.

Finite `SymbolicTensor::expand_dots` returns an `Atom` from the existing Spenso
component executor. It selects dots even inside opaque scalar payloads, retains
unsupported symbolic dimensions, and never realizes unrelated tensor factors.
Only intrinsic work is batched; sensitive callbacks retain their original
parse/execute and enclosing-normalizer order. This explicit component boundary
is separate from symbolic `undo_dots` and polynomial `expand`.

=== Trace materialization ownership

The R13 trace recurrence emits typed literal definitions through the existing
`SymbolicTensor` alias owner. Both ordinary and axial traces retain factored
subexpressions. `GammaSimplifySettings::expand_traces` and its builder are
removed; an explicit `expanded()` call uses the shared alias polynomial
materializer.

The previous trace-specific sparse emitter consumes an integer pairing recipe,
whereas the alias materializer accepts a general rational polynomial DAG,
including symbolic trace units and tensor-valued leaves. Reconstructing a
second recipe from that DAG would duplicate materialization and require another
eligibility dispatcher. The sparse emitter therefore remains only as a
test-side oracle, together with its overflow, allocation-bound and callback
regressions. Its 16 MiB row-buffer bound does not describe the shared alias
materializer. Exact comparisons between the oracle and explicit alias
materialization cover free and contracted words in four and symbolic
dimensions.

This resolves the sparse-emission entries below by consolidating their
production owner. It does not mark the remaining collection and contraction
entries complete; the final R17 audit still covers those paths separately.

== Direct production calls

#table(
  columns: (2.7fr, 1.5fr, 3.4fr),
  table.header([*Location at the fixed base*], [*Owner*], [*Replacement / disposition*]),
  [`idenso/src/dirac/simplify.rs:185`, `rewrite_node`],
  [Gamma rewrite, R13/R17],
  [Keep the rewritten trace body as an alias definition. Complete local identities without calling `expand`; surrounding factors stay opaque.],
  [`idenso/src/dirac/simplify.rs:839`, `evaluate_terminal_trace`],
  [Gamma rewrite, R13/R17],
  [Certify callback-produced bodies per distinct result and retain their sums in definitions. Never distribute them during callback cleanup.],
  [`idenso/src/dirac/simplify/trace_kernel.rs:624`, `evaluate_generic`],
  [Trace recurrence, R4/R17],
  [Return the shared recurrence's factored/aliased result. The sparse trace emitter remains only as an independent test oracle; explicit expansion uses the common alias materializer.],
  [`idenso/src/dirac/simplify/trace_kernel/sparse.rs:69`, `SparseOutput::emit`],
  [Trace materialization, R4/R17],
  [The bounded sparse recipe and its refusal behavior remain test-only coverage. Production explicit expansion uses the common alias materializer; its limits are independent of the test oracle row budget.],
  [`idenso/src/shorthands/schoonschip/contraction.rs:241`, local boundary product],
  [Contraction graph, R5/R12],
  [Visit one selected factor's alternatives, apply port overrides to neighbours, and merge identical remaining states. Delete the expanded local product.],
  [`idenso/src/shorthands/schoonschip/contraction.rs:742`, recursive network execution],
  [Contraction graph, R12],
  [Carry numerical weights on graph edges or alias definitions. Remove recursive scalar `expand_num` and the duplicate execution strategy.],
  [`idenso/src/shorthands/schoonschip/api.rs:246`, `NetworkSchoonschip::apply`],
  [Shared contraction orchestration, R12/R13],
  [Merge state coefficients directly; remove the post-pass `expand_num` and its expansion flag. R1 independently fixes the preceding sum accumulation.],
  [`idenso/src/lib.rs:647`, physical-leg adjoint],
  [Adjoint / gamma orchestration, R13/R17],
  [Attach boundary gamma-zero factors to each affected chain through the graph or an alias definition. Preserve the amplitude's scalar coefficients and factorization.],
  [`spenso/src/shadowing/collect.rs:129`, `Collectable::expand_with_map`],
  [Tensor collection, R12/R17],
  [Retain as caller-requested filtered expansion. Typed simplification uses the shared factor traversal instead; it does not invoke this explicit operation. Unselected coefficients remain opaque.],
  [`idenso/src/tensor/inference.rs:166,168,171,195`, `SymbolicTensor::expanded`],
  [Explicit materialization, R4/R14],
  [Retain variable-specific expansion and the typed-zero/interface policy. Route aliased polynomial expansion through the existing collector emitter; no simplifying verb invokes this method implicitly.],
  [`spenso/src/network/symbolica_interop.rs:288,355,385,426–427`, `TensorAtomMaps`],
  [Explicit component/network materialization],
  [The paired scalar/tensor `expand`, `expand_in`, `expand_num`, and `expand_via_poly` mappings remain caller-requested operations. They are not contraction scheduling.],
  [`spenso/src/tensors/parametric/atomcore.rs:1071,1080,1088,1099`, `TensorAtomMaps`],
  [Explicit component materialization],
  [Retain the component data mappings, with no new calls from symbolic contraction or inference.],
)

== Indirect distribution and materialization

Searching only for `.expand()` misses polynomial collection, the collector's
coefficient-list output, and trace emission. These are part of the same removal
boundary even though their spelling differs.

#table(
  columns: (2.7fr, 1.5fr, 3.4fr),
  table.header([*Location / operation*], [*Owner*], [*Replacement / disposition*]),
  [`idenso/src/selective_expand.rs`, `expand_in_patterns` and its metric / Minkowski / bispinor wrappers],
  [Contraction / explicit expansion, R12/R15],
  [Remove the pattern-search plus placeholder `coefficient_list` route. Shared graph incidence selects the affected component; explicit expansion alone emits its polynomial.],
  [`idenso/src/color/mod.rs:823`, `expand_color`],
  [Colour API, R13/R15],
  [Migrate callers to colour settings and explicit output materialization. Do not collect surrounding sectors merely to expose colour factors.],
  [`spenso/src/shadowing/collect.rs:180,191`, `collect_with_map`],
  [Graph-backed sector traversal, R12/R13],
  [Retain explicit filtered collection and original-slot frontend ingress. One shared helper returns the protected `AliasedAtom`; callers choose resolution. Typed domain passes use selected-factor traversal and their existing alias registry instead.],
  [`spenso/src/shadowing/collect.rs:247,253,478–505`, expansion dispatch],
  [Graph-backed sector traversal, R12/R15],
  [Delete the eight unused forwarding methods and migrate their callers. Keep the user-requested filtered operations that have a distinct explicit purpose; do not retain compatibility forwarding names.],
  [`idenso/src/shorthands/schoonschip/slot_contraction/components.rs:347–353`, coefficient-list emission],
  [Explicit aliased-result expansion, R5/R12],
  [Keep the tape/distributor and polynomial builder as the materialization owner. Contraction emits a typed alias DAG instead of invoking `to_expression` for every completed stage.],
  [`idenso/src/dirac/simplify/trace_kernel/sparse.rs:88–94`, coefficient-list emission],
  [Explicit trace materialization, R4/R17],
  [Retain integer coefficients and the existing checked row-buffer bound when the caller requests an expanded polynomial. Normal gamma simplification returns factored definitions.],
  [`idenso/src/dirac/simplify/trace_kernel.rs`, short-trace and axial output],
  [Gamma identities, R4/R13],
  [Identity-generated sums are valid definitions. Alias shared subresults and defer products-of-sums distribution, including the current implicit expanded axial terminal mode.],
  [`idenso/src/color/simplify.rs`, `collect_lines`, `collect_color`, `collect_rep_with_map`],
  [Colour orchestration, R6/R13],
  [Keep colour identities and projector substitutions. Run them on the shared graph, reusing unchanged branches and bulk construction; scalar and foreign tensor factors are opaque.],
)

== Additional audit sites

The post-M0 audit includes indirect distribution in addition to the direct
Symbolica calls above:

#table(
  columns: (2.5fr, 1.5fr, 3fr),
  table.header([*Operation*], [*Owner*], [*Required disposition*]),
  [`EXPANDSUMS` sum × sum helpers], [R12 contraction], [Delete with the network route; traverse one selected factor without constructing the distributed product.],
  [Default expanded short and axial 4D traces], [R4/R13 gamma], [Return alias recipes; retain expanded emission only at explicit materialization.],
  [Cofactor multiplication at collection sum roots], [R12 collection], [Keep the cofactor outside the root sum instead of copying it into each term.],
  [Colour payload distribution], [R6/R13 colour], [Preserve foreign factors and selected colour-sector aliases.],
  [Lazy tensor sums in component network execution], [Spenso network execution], [Keep lazy component sums; do not distribute symbolic graph numerators to schedule contraction.],
  [`SimplifySettings.expand`], [R14/R15 policy], [Remove implicit expansion from simplifying verbs.],
  [Generator grouping zero/equality test], [R16 production], [Use selected coefficient proofs; retain inconclusive results without expanding the graph numerator.],
  [Factorized `factor_terms`], [R5 contraction], [Stream alternatives while keeping branch-local dummy scope; copied contracted subexpressions cannot share dummy identities.],
  [Alias materialization], [R4 explicit materialization], [Use local polynomial variables and bulk construction; keep the already-expanded shortcut. This operation is never called implicitly by contraction.],
)

== Historical policy switches

At the fixed base, `GammaSimplifySettings::expand_traces` controlled trace-body cleanup,
terminal admission, output choice and sparse dispatch. Removing only its final
`expand` call would leave those policies inconsistent. R4/R13 replace its
output role with explicit aliased-result materialization; the callback and
closure checks remain with inference. R15 migrates Rust and Python callers,
the `with_expanded_traces` builder, stubs, tests and notebook examples together.

`SchoonschipSettings::expand_contracted_sums` and the network route's
`EXPANDSUMS` const parameter influenced local contraction, collector
admission, recursive strategy selection and scalar coefficient expansion.
R5/R12 replace their contraction role with one-factor traversal. They are
removed with the old strategies, and R15 updates the surface in one change.

== Preservation and verification

The new engine must preserve unchanged atoms, logical port order, unresolved
port occurrences, typed zeros, metadata, disconnected components and scalar
spectators. Callback-created results still need certification: the existing
`g(a,b) * T(a)` normalizer that turns `T(b)` into a scalar prevents blindly
retaining the predicted interface. Alias handles with relabeled ports need
one registered definition per distinct labelling until Symbolica supports
parametric aliases. None of these checks is replaced by expansion.

Search scope is all Rust sources under `crates/idenso/src` and
`crates/spenso/src`, including `expand_in`, selective wrappers,
`collect_symbol`, `coefficient_list`, and polynomial emission. Test-only
modules and the `TensorAtomMaps` trait declarations were inspected separately;
they do not add production call sites. In particular,
`reference_cases.rs:485` and `tensors/parametric.rs:3330` are test assertions.

Each removal is checked by the R3 harness and existing exact, callback and
HEP-component suites. Synthetic materialization fixtures can compare explicit
expansion with FORM; captured graph numerators use factor-preserving or exact
component checks, with the latter's sampled momentum and dimension scope
recorded explicitly. The final audit repeats
this inventory and checks that the remaining calls are explicit materialization
or tests, rather than hiding distribution under a renamed helper.

== Audit additions at `15d64dd1` and R5

The following indirect sites were missing from the initial inventory. They are
requirements for the same R17 removal gate, not an additional expansion API.

#table(
  columns: (2.5fr, 1fr, 3.7fr),
  table.header([*Owner / site*], [*Register*], [*Required disposition*]),
  [`EXPANDSUMS` sum-by-sum helpers in the legacy network contractor], [R12],
  [Delete with that engine after migrating its UV caller. No distributed-Atom fallback in the shared contractor.],
  [`dirac/simplify.rs`, default short 4D and axial trace terminals], [R13],
  [Emit local identities into the existing alias registry, preserving repeated-index scope. Sparse expanded output is a test oracle only.],
  [`shadowing/collect.rs`, root-Add cofactor copying], [R12/R13],
  [Keep the cofactor outside selected alternatives. Carry selected sectors through the shared tape and alias owner.],
  [Colour payload distribution and local terminal sums], [R13],
  [Visit selected alternatives without multiplying foreign or scalar spectators into each term.],
  [`network/contract.rs`, lazy tensor sums], [R8/R16],
  [Keep symbolic alternatives in the structural graph until explicit component execution. Preserve the component execution boundary.],
  [`SimplifySettings.expand`], [R13/R15],
  [Remove; only the explicit materialization verb may request expansion.],
  [`feynkit-generator/grouping.rs`, colour-group zero tests], [R16],
  [Check each colour coefficient separately; do not distribute the complete graph numerator.],
  [Factor tape powers and in-flight `factor_terms`], [R5/R12],
  [Retain alternatives and scope internal dummy indices per copy before contracting surviving ports. Avoid enumerating ordered equivalent tuples.],
  [`aliases/materialization.rs`, polynomial forward pass], [R4/R5],
  [Explicit output boundary only. Keep the expanded shortcut, use definition-local variables and one bulk additive merge.],
)

The power regression oracle contracts the base before FORM powers or evaluates
its components directly. Expanding `(p(a)*q(a)+x)^2` first destroys the scope
information under test and is not an independent oracle.

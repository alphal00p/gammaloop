= Tensor engine expansion inventory
<tensor-engine-expansion-inventory>

This is the M0 inventory for R17 of
#link("tensor-engine-consolidation-plan.typ")[the tensor engine consolidation plan].
It records the fixed implementation base `bfb80d8c`, including the preserved
ladder benchmark and bounded sixteen-gamma trace work. It is an inventory of
work to do, not a claim that the internal expansion paths have been removed.
Line numbers below refer to that base; function names identify the owners
when subsequent changes move them.

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
responsibility. Expanding a difference in a correctness assertion is also
outside the production engine and remains a useful exact check.

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
  [Return the shared recurrence's factored/aliased result. The sparse emitter remains available only to explicit materialization.],
  [`idenso/src/dirac/simplify/trace_kernel/sparse.rs:69`, `SparseOutput::emit`],
  [Trace materialization, R4/R17],
  [A row-budget refusal retains the completed recipe; explicit expansion may walk its terms. Do not replace a refused allocation with an internal Atom expansion.],
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
  [Replace `expand_in(COLLECT)` with the shared factor traversal. Selected sectors are visited locally; unselected coefficients remain opaque.],
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
  [Replace `collect_symbol(COLLECT)` distribution with selected-factor iteration and shared alias definitions. Preserve its useful protection of scalar spectators.],
  [`spenso/src/shadowing/collect.rs:247,253,478–505`, expansion dispatch],
  [Graph-backed sector traversal, R12/R15],
  [Migrate `expand_chains_and_traces`, `expand_tensors`, `expand_tagged_tensors`, `expand_rep`, `expand_reps`, `expand_metrics`, and their mapping callers to graph traversal. Do not retain forwarding compatibility names.],
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

== Policy switches to remove

`GammaSimplifySettings::expand_traces` currently controls trace-body cleanup,
terminal admission, output choice and sparse dispatch. Removing only its final
`expand` call would leave those policies inconsistent. R4/R13 replace its
output role with explicit aliased-result materialization; the callback and
closure checks remain with inference. R15 migrates Rust and Python callers,
the `with_expanded_traces` builder, stubs, tests and notebook examples together.

`SchoonschipSettings::expand_contracted_sums` and the network route's
`EXPANDSUMS` const parameter currently influence local contraction, collector
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
HEP-component suites. Its factorized result is compared only after an explicit
test-side expansion with the FORM-certified polynomial. The final audit repeats
this inventory and checks that the remaining calls are explicit materialization
or tests, rather than hiding distribution under a renamed helper.

= Tensor engine consolidation plan
<tensor-engine-consolidation-plan>
#quote(block: true)[
#strong[Lifecycle:] Design proposal and plan of record. No code in this plan
has been written; the measurements are point-in-time evidence taken on
2026-09-27 against `a27ca81d8958` with release wheels of the Python
bindings, and they must be re-taken by the benchmark harness in M0 before
they are used as acceptance thresholds.

#strong[Scope:] the index-contraction and simplification stack shared by
`spenso`, `idenso` and `spynso3` (`TensorExpression`), and the way
`gammalooprs` consumes it. Concrete tensor-component execution stays with
Spenso's tensor libraries and is touched only where it must consume the new
result type.

#strong[Inputs:] the API audit of `TensorExpression` and its production
callers, the gluon-ladder and gamma-trace comparisons with FORM, the
keep-the-rest probes, the partial-parse and network cost measurements, the
factor-graph contraction prototype, and the internal-path audit of the
second agent working on this stack (reproduced in the register below).
]

== Decision in one paragraph
<decision>

There will be one contraction engine, and it will run on the network. The
network's structural layer becomes cheap (borrowed leaves, owned incidence,
port overrides) and total (every factor it does not own is an opaque leaf that
passes through unchanged). Feynman-rule application, contraction and
materialization become three separately timed steps: rules are applied by
Symbolica on the factorized expression, contraction expands one factor at a
time and pushes its index structure into the still-factorized neighbours
while merging states whose remaining factors are identical, and the result is
an aliased tensor expression, a DAG whose nodes are typed tensor aliases,
expanded to a polynomial or a plain atom only on demand. Rules and tensors are
certified once and carry their certificates; nothing is re-proved on the
production path. The public surface shrinks to a small set of verbs, and the
internal duplicates identified by both audits are retired under the existing
regression suites.

== Evidence
<evidence>

All timings: release wheels, medians of alternating rounds where noted, one
host; FORM 5.0.0 as the external reference. "Module" is the fused
replace-and-contract route ported in the `tensor-module` change; it is a
stepping stone and is superseded by this plan.

=== Gluon ladder (8 three-gluon vertices, 11 dummies, scalar result of 9,652 terms)

#table(
  columns: (2.6fr, 1fr, 1fr, 2.2fr),
  table.header([*Route*], [*Original order*], [*Reverse order*], [*What it measures*]),
  [typed accumulator (notebook)], [2.71 s], [—], [replace, multiply, expand, contract per vertex],
  [`replace_tensor` + `expand` + `schoonschip`], [0.94 s], [0.35 s], [materialize generated terms, re-parse],
  [`replace_tensor` + collector (`expand_contracted_sums`)], [0.84 s], [0.29 s], [current fast route],
  [scalar Symbolica recipe (control)], [0.86 s], [0.32 s], [pure Symbolica replacement rules],
  [module `replace_tensor(contract=True)`], [0.54 s], [0.21 s], [no materialized intermediate],
  [FORM], [0.73 s], [0.17 s], [reference],
  [factor-graph DP, contraction only, Python prototype], [1.05 s], [*0.107 s*], [1,730 edge expansions in reverse order],
  [factor-graph DP, Feynman rules for all vertices at once], [0.4 ms], [0.4 ms], [`replace_multiple`, factorized],
  [factor-graph DP, polynomial forward pass to the expanded result], [—], [0.10 s], [only when the polynomial is wanted],
)

Generated versus retained terms in the expanded routes: 50,684 / 10,947 at
stage 6, 92,341 / 23,937 at stage 7, 186,516 / 9,652 at stage 8. The DP
performs 17,515 edge expansions in the original order and 1,730 in the
reverse order, with 1,320 and 105 peak live states, and its result is a
59 KB DAG against a 1.9 MB polynomial. All routes agree with the
FORM-certified polynomial.

=== Gamma traces

The Dirac algebra is at parity with or ahead of FORM (4D free 12 and 14
gammas: 3.1 and 25 ms against `trace4` 3.5 and 21 ms; D-dimensional: 2.1 and
45 ms against `tracen` 5.1 and 65 ms). FORM's expanded output costs six to
fifteen times more only because `expand()` materializes every term at four to
seven microseconds each. Streaming that expansion through the collector is
correct and not faster; there is nothing to merge.

=== Cost of discovering that there is nothing to do

#table(
  columns: (2.6fr, 1fr, 1fr, 1fr),
  table.header([*Idempotent call on a contracted stage output*], [*6,069 terms*], [*23,937 terms*], [*9,652 scalar terms*]),
  [`schoonschip()`], [64 ms], [124 ms], [15 ms],
  [`simplify_metrics()`], [35 ms], [74 ms], [15 ms],
  [`TensorExpression(expr)` from a bare atom], [205 ms], [231 ms], [116 ms],
  [module admission and traversal, no match], [6 ms], [22 ms], [9 ms],
)

A no-op pass costs as much as the productive pass that produced its input,
because every call re-classifies every leaf. Production works on bare atoms
(`simplify_metrics` 25 call sites, `to_dots` 28, `normalize_dots` 19,
`simplify_gamma` 13, `simplify_color` 15, `replace_multiple` 171,
`schoonschip` 2, `schoonschip_net` and `replace_tensor` 0), so this is paid on
every call. `simplify_metrics` also disables rank-one handling, which bypasses
every collector fast path.

=== Network versus collector, per ladder-shaped term (Rust, release)

#table(
  columns: (3fr, 1fr, 1fr),
  table.header([*Step*], [*idempotent*], [*productive*]),
  [collector `schoonschip()`, whole pass], [5.7 µs], [6.3 µs],
  [partial symbolic parse (`depth_limit 1`, opaque fast inference)], [23.3 µs], [28.4 µs],
  [symbolic parse without depth limit], [65.8 µs], [71.4 µs],
  [`schoonschip_net()` = parse, execute, normalize], [169 µs], [186 µs],
  [`canonize` (graph canonization)], [430–520 µs], [—],
)

The parse is 15% of the network route and already four times the collector's
whole pass; execution is 85%, spent multiplying opaque leaves into one result
tensor even when nothing contracts. The whole-sum figures previously quoted
for `schoonschip_net` (0.7–0.9 ms per term) were dominated by an O(n²)
accumulation in its driver, not by parsing.

=== Keep-the-rest conformance

`simplify_color`, `simplify_gamma`, `simplify_metrics`, `schoonschip` and the
collector route keep `(1+x)·(A + f(z)·B)` factored and pass foreign factors
through; the current fast route even keeps disconnected components factored.
The module distributes scalar spectators and refuses any foreign factor. One
defect: a metric contracting into a gamma with `AUTO` spinor indices fails in
`simplify_metrics` and `simplify_gamma` while explicit indices work.

== Work in flight on the second agent's line
<second-agent>

Landed on the shared base (`e116e3c8` and its parents): reuse of contraction
state across Schoonschip normalization stages, preserved diagram indices and
scoped tensor indices, tensor-library membership resolution, edge
denominators and shared propagator metadata, and MiTeX rendering of model
parameter names. In flight and uncommitted at the time of writing (about
4,700 lines across `examples/notebooks`, `examples/reproducers` and
`docs/products/idenso`): the multi-loop FORM comparison — a massless
three-loop propagator with a fermionic outer ring, then four-loop fermionic
and gluonic ladders with physical momentum routing and a metric projection of
the external indices — with a shared benchmark driver (`fermion_ladder.py`,
`gluon_ladder.py`), native FORM programs, staged timing records and an exact
check battery (polynomial equality with FORM, D → 4 specialization, 72 exact
HEP-component comparisons).

Their routes matter to this plan. The gluonic case already works outside-in:
tensor-safe vertex replacement, local Schoonschip contraction and contraction
against the remaining opaque vertices in a measured order (`1,2,8,3,7,4,6,5`
against FORM's `5,4,6,3,7,2,8,1`), without expanding the product of eight
vertices first. The fermionic route runs the typed gamma simplifier to a
scalar trace, converts once to a polynomial, substitutes momenta group by
group with a collect after each, and emits one Atom at the end; against the
previous simultaneous-substitution driver this gave paired CPU reductions of
27.5% (three-loop D), 18.2% (four-loop 4D) and 53.9% (four-loop D).

#table(
  columns: (1.5fr, 0.6fr, 1fr, 0.9fr, 0.7fr, 1fr, 0.9fr, 0.7fr),
  table.header([*Numerator*], [*Dim.*], [*Idenso, earlier*], [*FORM*], [*Ratio*], [*Idenso, current*], [*FORM*], [*Ratio*]),
  [Three-loop fermion], [4D], [1.027], [0.505], [2.03×], [0.612], [0.356], [1.72×],
  [Three-loop fermion], [D], [4.414], [2.560], [1.72×], [2.879], [2.231], [1.29×],
  [Four-loop fermion], [4D], [13.032], [8.875], [1.47×], [6.523], [6.583], [0.99×],
  [Four-loop fermion], [D], [231.908], [60.250], [3.85×], [133.257], [71.500], [1.86×],
  [Four-loop gluon], [4D], [664.375], [622.000], [1.07×], [773.030], [778.000], [0.99×],
  [Four-loop gluon], [D], [1,068.251], [1,089.000], [0.98×], [1,410.728], [1,479.000], [0.95×],
)

Median CPU milliseconds of the complete algebra (trace or vertex replacement,
contraction, momentum routing, final scalar collection), five alternating
warm batches on one pinned CPU, FORM 5.0.0; ratios divide medians. "Earlier"
is the cohort reported during the working session, "current" the cohort in
the notebook text dated 2026-09-27 after the routing change; the gluonic
algorithm did not change between them, so its rows show the host's drift.

The remaining gap is the four-loop D-dimensional fermion trace. Their
single-call diagnostic splits it into 88.6 ms tracing (including collection
of the trace into scalar form), 8.8 ms converting to a polynomial, 29.7 ms
substituting momenta and 7.4 ms emitting the final scalar.

What this plan takes from their work:

- The M0 harness is their driver extended with the ladder, trace and
  keep-the-rest cases, not a new one; the four-loop D fermion case is its
  headline, and their exact check battery is the correctness gate of every
  milestone.
- The 46 ms after the trace are materialization steps this plan makes
  on-demand: with an aliased result the routing substitutions act on alias
  definitions and the polynomial is built once, if at all.
- The tracing step includes collecting the trace into expanded scalar form;
  the Dirac algebra itself measured at parity with `tracen` in this session.
  The aliased result type lets the trace stay factored into the routing step.
- Their outside-in gluonic route and the DP contraction differ only in
  merging states by remaining forms; R5 is measured against their gluonic
  numbers as well as the historical ladder.

== Principles
<principles>

+ #strong[Keep the rest.] A verb acts on the structure it owns and returns
  everything else exactly as found: scalar spectators and factored sums stay
  factored, foreign tensors pass through, disconnected components are never
  multiplied out. Declining is acceptable for malformed input only, never for
  unfamiliar input. Colour and gamma simplification never require a full
  expansion.
+ #strong[Three steps, three clocks, and no internal expansion.] Rule
  application (Symbolica, factorized), contraction (the engine) and
  materialization (polynomial, atom, evaluator) are separate operations with
  separate benchmarks. `expand()` is called only by the user, for example to
  compare with FORM; no verb, option or internal step expands an expression,
  neither as its output form nor as an intermediate (decision D5, item R17).
+ #strong[Own once.] Source atoms are borrowed for the life of an operation;
  rule outputs are owned once per distinct binding; the graph is owned and
  arena-reused; atoms are created at emission, one per distinct variable, and
  never per generated term.
+ #strong[Trust once.] A rule is compiled and certified at construction. A
  tensor produced by the engine carries its interface and a normal-form
  marker, so admission and idempotent calls are O(1). Bare atoms pay inference
  once, at construction. The hard constraint from the second audit stands: a
  callback right-hand side can change the interface a contraction predicted
  (`g(a,b)·T(a)` where normalization turns `T(b)` into a scalar), so callback
  outputs are certified per distinct output, and only pattern right-hand sides
  are certified structurally once per rule.
+ #strong[One graph.] The partial parser's structure is the only graph.
  Contraction, interface inference, canonization, component partition,
  dummy allocation and dot conversion read it; nothing re-walks atoms to
  rediscover it.
+ #strong[Results are aliased tensor expressions.] The engine returns a root
  plus typed alias definitions. Expansion is a method, not a side effect.
+ #strong[No new frameworks.] Extend `SymbolicTensor`, `SlotMatcher`,
  `InterfaceInference`, Spenso's `Network` and Symbolica's `AliasedAtom`;
  do not add another tensor type or a generic transformation layer.

== Target architecture
<target-architecture>

=== Symbolica

- Matching, conditions, wildcard binding and right-hand-side instantiation are
  Symbolica's and remain so. `replace` as an expression rewriter is used for
  rule application on the factorized expression; it is not used to drive
  contraction, because its output is a tree that must be expanded and
  re-parsed.
- `AliasedAtom` (root, alias table, nested definitions, `register_alias`,
  fused arithmetic, `into_inner`, `evaluator_multiple` via
  `EvaluatorBuilder::add_aliases`) is the result container. It needs one
  extension: parametric aliases, where a handle with slot arguments resolves
  by pattern rather than by literal atom equality, so a tensor alias used with
  relabeled or vector-filled ports resolves to its definition with the
  corresponding substitution. Until that lands, the engine registers one
  literal alias per distinct use-labelling.
- Requested upstream, not blocking: a streaming `replace`/`expand` that hands
  generated terms to a consumer without assembling the sum.

=== Spenso: the structural network layer

- A borrowed structural graph: leaves are `AtomView`s (into the source or into
  the owned right-hand-side arena) with a port list from the fast syntactic
  inference and an incidence table; the graph is arena-allocated and reused
  across terms and operations. Parse settings are the existing partial ones
  (`depth_limit`, `ShorthandParsing::Opaque { Fast }`, pre-contracted
  scalars); the strict tensor filter decides which heads are leaves.
- Two kinds of leaves. Library leaves (metric, identity, tagged vectors) have
  algebra; every other leaf, including gammas and colour structures in a
  metric pass and alias handles everywhere, is opaque. Contraction rewrites
  only edges with a library endpoint. An opaque leaf's ports may be relabeled
  or receive a compact vector (a port override); its content is never read and
  a new atom for it is built once, at emission.
- Execution semantics are outside-in and local: replacing a leaf splices the
  right-hand-side template port by port and contracts only the edges the
  splice exposed. There is no reduction to a single result tensor and no
  global scheduling; the caller's replacement order plus the contraction graph
  decide the order.
- Materialization is the existing owned `Network` and `NetworkStore`, built
  from the structural graph when component data are needed. The scalar-store
  aliases already produced there become the same alias handles the engine
  emits.
- Returned inference results keep logical port ordering (second audit, item
  2), and dummy names are reserved once per operation.

=== Idenso: the engine

- #strong[Rule application.] All rules of a stage at once, on the factorized
  expression, through certified `TensorRule` objects (pattern compiled once,
  wildcard closure checked once, pattern right-hand sides certified
  structurally once, callback right-hand sides certified per distinct output
  and cached by binding).
- #strong[Contraction.] Read the contraction graph off the factors; choose an
  order from it; expand one factor at a time; push each term's metric and
  vector structure into the neighbours by port override; resolve a contraction
  that becomes internal to one factor with the per-term reducer; merge states
  whose remaining factors are identical, adding their coefficients. The
  per-term reducer is the collector's `reduce_components` logic on the
  structural graph; the collector's tape, distributor and coefficient-list
  emitter remain the expansion-on-demand path.
- #strong[Emission.] One alias per state, one definition per alias, a root
  handle; opaque leaves with port overrides materialized once per distinct
  variable.
- #strong[Other passes.] Gamma, colour and epsilon keep their algebra and run
  on the same graph; they alias what they do not own instead of expanding
  around it; shared orchestration carries observations between passes and
  invalidates them on change (second audit, item 4).

=== Spynso3: the surface

- `TensorExpression` stays the typed atom. `AliasedTensorExpression` is added:
  root `TensorExpression`, aliases as `(handle, body)` pairs of
  `TensorExpression`s whose interfaces match, with `to_expression()`
  (resolve), `expand()` (polynomial forward pass), `evaluator(...)`,
  `map_aliases(f)`, and every verb acting on the root with handles opaque.
  The rank-zero case is the scalar `AliasedExpression`; no separate wrapper.
- `TensorRule` is exposed so a rule is built once and applied many times.
- Transformation policy (identity results, typed zeros, inference, lost
  interfaces) has one Rust owner in `SymbolicTensor`; the bindings convert,
  describe, wrap and translate errors (second audit, item 6).

== Public verbs
<public-verbs>

Target surface of `TensorExpression` for index algebra, with the current names
each one absorbs. Per decision D1 there is no compatibility layer: the renames
land in one change together with every caller.

#table(
  columns: (1.2fr, 2.6fr, 2.6fr),
  table.header([*Verb*], [*Meaning*], [*Absorbs*]),
  [`replace(rule | rules)`], [tensor-safe rule application, factorized output, certified once per rule], [`replace_tensor`, `replace_multiple` for tensor leaves, `contract=` flag],
  [`contract(order=None, output="aliased")`], [metric, identity and vector contraction on the structural graph; output aliased, or `"expanded"` on request], [`schoonschip`, `schoonschip_net`, `simplify_metrics`, `expand_metrics`, `expand_mink`, `expand_bis`, `expand_mink_bis`, `expand_in_patterns`, `undo_schoonschip`],
  [`simplify_gamma(settings)`], [unchanged algebra, aliases foreign structure], [`simplify_gamma0`, `simplify_gamma_conjugate` become settings],
  [`simplify_color(settings)`], [unchanged algebra, aliases foreign structure], [`expand_color`, `collect_color`, `collect_color_constants`, `wrap_color` become outputs or settings],
  [`simplify_epsilon()`], [unchanged], [—],
  [`canonize()`], [graph canonization on the shared graph], [`wrap_dummies`, `wrap_indices` where they only served canonization],
  [`to_dots()` / `undo_dots()`], [representation toggles; distinct paths by design], [`normalize_dots`, `metric_shorthand_to_dot`, `expand_dots`],
  [`expand()` / `evaluator()`], [materialization, explicit], [—],
)

Retained as distinct because they handle different situations (second
audit): intrinsic `g` normalization, nested-vector normalization and product
contraction; the bulk collector and the callback-sensitive fallback; sparse
and factored trace outputs; a network's semantic expression, materialized
indices and component data. Verbs with no caller outside tests and the
notebook (about a quarter of the current ninety-odd methods) are removed in
the same change unless a product doc claims them.

== Consolidation register
<register>

Items are grouped by the milestone that delivers them. "Second audit" marks
findings of the other agent; their line references are to
`crates/spynso3/src/network.rs` (`multiply_network`),
`crates/idenso/src/tensor/inference.rs` (fallback inference),
`crates/idenso/src/tensor/replacement.rs` (`Signature::observe`),
`crates/spynso3/src/expression.rs` (`simplify` pipeline,
`from_transformed_atom`), `crates/idenso/src/color/simplify.rs`
(`rewrite_terms`, `ProductView`) and
`crates/idenso/src/shorthands/schoonschip/normalize_dots.rs` (preflight).

#table(
  columns: (0.5fr, 3fr, 1.6fr, 0.8fr),
  table.header([*Id*], [*Item*], [*Source*], [*Milestone*]),
  [R1], [O(n²) `sum +=` accumulation in the network schoonschip driver], [measured], [M0],
  [R2], [metric into gamma with `AUTO` spinor indices fails], [measured], [M0],
  [R3], [benchmark harness: ladder, traces, keep-the-rest probe, idempotent cost, production-derived numerator; alternating interpreters, medians, FORM reference], [this plan], [M0],
  [R4], [`AliasedTensorExpression` and `TensorRule` types; parametric aliases in Symbolica or literal-per-labelling fallback], [this plan], [M1],
  [R5], [factor-graph DP contraction in Rust over the existing collector graph, `contract(output="aliased")`], [prototype], [M1],
  [R6], [colour `rewrite_terms` and `ProductView` rebuild sums and products by repeated binary operations; colour retries a failed root rewrite], [second audit, 5], [M1],
  [R7], [dot-normalization preflight constructs a replacement the rewrite constructs again], [second audit, 5], [M1],
  [R8], [borrowed structural layer in Spenso: leaf views, ports, incidence, port overrides, component partition, arena reuse; materialization to the owned network], [this plan], [M2],
  [R9], [Spenso returns logical port ordering; Idenso stops rescanning arguments; dummy names reserved once], [second audit, 2], [M2],
  [R10], [one product plan and one batched relabel for composition; logical composition in `SymbolicTensor`, graph and store maintenance in Spenso], [second audit, 1], [M2],
  [R11], [`Signature::observe` and scope validation share one traversal through `InterfaceInference`; certified rules and carried certificates make the per-call source proof disappear], [second audit, 3; this plan], [M3],
  [R12], [engine on the structural layer: rules at once, DP contraction, per-term reducer, emission; retire `contraction.rs` strategies and the graph half of `slot_contraction`/`components`], [this plan], [M3],
  [R13], [shared Rust orchestration of passes with carried observations and invalidation; Python pipeline stops nesting cleanup loops; gamma and Schoonschip stop rescanning unchanged input], [second audit, 4], [M3],
  [R14], [one owner for transformation policy (`from_transformed_atom` duplicates `SymbolicTensor`'s decisions)], [second audit, 6], [M3],
  [R15], [public verbs: renames and removals in one change with every caller, `.pyi` and product docs regenerated; no compatibility aliases], [audit; decision D1], [M4],
  [R16], [`gammalooprs` moves from `simplify_metrics` on bare atoms to certified `contract` with aliased output; evaluator consumes aliased results; UV profile re-measured], [audit], [M5],
  [R17], [no internal expansion: inventory every production `expand`, `expand_num`, `expand_via_poly` and `expand_in_patterns` call in `idenso` and `spenso` (network driver, colour coefficient collection, trace emission, collector output) and replace each with term iteration that builds no atoms; the options `expand_traces` and `expand_contracted_sums` go], [decision D5], [M0 inventory; M1–M3 removal],
  [R18], [standalone `rust-script` and Python reproducers under `examples/reproducers/symbolica-aliases/` showing why each Symbolica request is needed: parametric aliases, streaming replace/expand, evaluators from aliased atoms in Python], [decision D4], [M0],
)

== Milestones and gates
<milestones>

Every milestone lands under the existing exact-identity and HEP-library
regression suites plus the M0 harness; a milestone is accepted when its gates
hold on the harness, not on ad-hoc runs. Ownership: the second agent's line is
the base; parser and collector files under their edit land first, and R8 and
R12 are designed with them before code is written. Workflow (decision D3):
every register item is its own jj change created with `jj new default@-` when
the item starts, so its parent is the second agent's last completed change and
does not move while the item is in progress.

+ #strong[M0 — measure and stop the bleeding.] R1, R2, R3, R17 (inventory),
  R18. Gates: the network driver's whole-sum per-term cost within 2× of its
  single-term cost; the `AUTO` reproducer passes; the harness reproduces the
  evidence tables above within host drift; the expansion inventory names an
  owner and a replacement for every call site; each Symbolica reproducer runs
  and prints the limitation it demonstrates.
+ #strong[M1 — the result type and the algorithm, on today's graph.] R4,
  R5, R6, R7. Gates: `contract(output="aliased")` on the ladder equals the
  FORM-certified polynomial after `expand()`; contraction-only time in the
  reverse order under 50 ms and the aliased result under 100 KB; the
  polynomial forward pass under 150 ms; colour and dot changes are
  behaviour-preserving on their suites.
+ #strong[M2 — the structural layer.] R8, R9, R10. Gates: partial structural
  parse per ladder-shaped term within 2× of the collector's whole pass;
  interface inference and `canonize` read the graph without a second walk;
  composition builds one plan.
+ #strong[M3 — one engine.] R11, R12, R13, R14. Gates: `contract` on an
  already-contracted certified input is O(1); no verb calls another verb's
  full pass on unchanged input; the slot engine's strategies and graph code are
  deleted; ladder and trace numbers do not regress against M1; every
  keep-the-rest probe passes, including foreign factors and factored sums
  inside the fused replace-and-contract path.
+ #strong[M4 — the surface.] R15. Gates: the public-API test lists the
  target verbs and nothing else; every caller in production, notebooks and
  tests uses them; docs and stubs regenerate cleanly.
+ #strong[M5 — production.] R16. Gates: `gammalooprs` numerator pipelines
  use certified rules and aliased results; the aa→aa three-loop smoke and the
  UV scalar profile are re-measured and recorded.

The `tensor-module` change is not merged; its tests for fused
replace-and-contract carry over to `replace(...).contract()` in M1, and its
sparse-emission fallback carries over to the expansion path.

== Decisions
<decisions>

Taken on 2026-09-27 unless marked open.

- #strong[Direction and order.] Accepted with the answers below: one engine
  on the network, aliased tensor expressions as results, the DP contraction,
  M0 to M5 in order, the second agent's line as the base.
- #strong[D1 — old Python names: no compatibility layer.] The renames land in
  one change together with every caller: the production sites
  (`simplify_metrics` 25, `to_dots` 28, `normalize_dots` 19, `simplify_gamma`
  13, `simplify_color` 15, and the tensor-leaf uses among 171
  `replace_multiple`), the notebooks, the tests, and the generated stubs and
  docs. Nothing is kept for compatibility.
- #strong[D2 — `TensorNetwork` in Python: open.] Facts for the call: no
  production Python caller (`gammaloop-api` and `feynkit-py`: none); 12
  notebook and 13 product-doc mentions and 26 test uses, all for execution and
  display (`to_network`, `execute`, `step`, `result_tensor`, `to_tensor`,
  `result_scalar`, `to_dot`, `render`, `to_linnest`); `gammalooprs` uses
  networks in Rust in six files, for execution. Recommendation: keep it as the
  execution and rendering object and remove its simplification entry point
  (`schoonschip_net`); removing the type would only move those same methods
  onto `TensorExpression`.
- #strong[D3 — separate changes on a fixed base.] Every register item is its
  own jj change created with `jj new default@-` when the item starts: its
  parent is the second agent's last completed change, not their working copy,
  so the base does not move while the item is in progress. This plan change is
  rebased onto that base.
- #strong[D4 — Symbolica requests come with reproducers.] Each request is
  accompanied by a standalone `rust-script` or Python minimal example under
  `examples/reproducers/symbolica-aliases/` that shows the limitation and why
  the engine needs the change: parametric aliases, streaming replace/expand,
  evaluators from aliased atoms in Python (R18).
- #strong[D5 — expansion is the user's, never the engine's.] `expand()` is
  called only by the user, for example to compare with FORM. No verb, option or
  internal step expands an expression, neither as its output form nor as an
  intermediate; R17 removes the existing internal expansions.

== Risks
<risks>

- The structural layer's constant factors: the collector's 6 µs per term is a
  flat incidence table with union-find; a generic graph may not reach it.
  Mitigation: M2's gate is 2×, and the per-term reducer keeps its flat form.
- State growth in the DP with a poor order (17,515 edges against 1,730 on the
  ladder). Mitigation: order chosen from the contraction graph; budgets that
  flush partial results instead of failing.
- Callback right-hand sides can change predicted interfaces. Mitigation:
  certification per distinct callback output remains; only pattern rules are
  certified once.
- Alias resolution with port overrides is only correct with parametric
  aliases or literal-per-labelling registration. Mitigation: the latter is the
  M1 fallback and is tested by resolving and comparing with the expanded
  route.
- Two agents editing the same files. Mitigation: D3, and R8/R12 designed
  jointly before code.
- Removing every internal expansion (R17) may surface algorithms that only
  work on expanded input today (trace emission, colour coefficient
  collection, the collector's polynomial output). Mitigation: the M0
  inventory sizes each site before M1 commits to a replacement; term
  iteration over factorized input is the default replacement.

== Measurement protocol
<protocol>

Alternating interpreters per round, medians of three rounds, a pure-Symbolica
scalar recipe as the control column, FORM as the external reference, release
wheels only, and per-stage term or edge counts recorded next to times. The
prototype scripts (`factor_dp.py`, `factor_dp_dag.py`,
`factor_dp_poly.py`, `validation_overhead.py`, `idempotent_cost.py`,
`partial_parse_cost.py`, `keep_the_rest_probe*.py`) move from the session
scratchpad into the repository as the M0 harness.

== Related documents
<related>

- #link("schoonschip-net-parsing.typ")[Schoonschip network architecture] —
  the current network path this plan replaces.
- #link("network-simplification-status.typ")[Network simplification status] —
  earlier measurements of the same path.
- #link("idenso-architecture.typ")[Idenso implementation architecture] and
  #link("spenso-architecture.typ")[Spenso architecture] — the layers this plan
  changes.
- #link("api-documentation-debt-register.typ")[API documentation debt register]
  — documentation obligations for the renamed surface.

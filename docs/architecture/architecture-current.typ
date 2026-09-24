= GammaLoop architecture
<gammaloop-current-architecture>
#quote(block: true)[
#strong[Reviewed:] 2026-08-18 against the verified source scopes registered in
`docs/developers.toml`.

#strong[Lifecycle:] Current implementation architecture. The review
checked the state, evaluator, event-processing, integration, and
quantity contracts named below against their implementations.
]

== Scope
<scope>
This document describes the current, implemented architecture of the
GammaLoop application. It is not an architecture overview of the entire
multi-product repository; Linnet, Spenso, Idenso, and Vakint have their
own manuals and reference surfaces.

The GammaLoop application is organized around two principal Rust crates
while using those shared workspace libraries:

- `gammalooprs` (`crates/gammalooprs`): core physics/domain logic, graph
  processing, integrand construction, evaluation, and integration.
- `gammaloop-api` (`crates/gammaloop-api`): CLI, REPL, Python bindings,
  command parsing, and persisted state orchestration.

== High-Level Architecture
<high-level-architecture>
At a high level, gammaLoop uses a layered architecture with a stateful
application shell.

```text
User / Automation
  -> CLI binary (`gammaloop`) and REPL
  -> Python module (`_gammaloop` via pyo3)
  -> run cards (`run.toml`)

Application Layer (`gammaloop-api`)
  -> clap command parsing
  -> command execution against mutable `State`
  -> state/load/save and logging control

Domain Layer (`gammalooprs`)
  -> model + parameter cards
  -> graph and process representations (amplitude/cross-section)
  -> preprocessing and CFF generation
  -> GL integrand construction + evaluator stacks
  -> differential event generation, selectors, and observables
  -> Monte Carlo integration and stability checks

Infrastructure
  -> filesystem persistence (`gammaloop_state/`)
  -> optional local scratch persistence (`.local/scratch/...`)
  -> Symbolica state import/export
  -> tracing/logging
```

== Main Components
<main-components>
=== 1. Entry Points and Interfaces
<1-entry-points-and-interfaces>
- CLI entry point: `crates/gammaloop-api/src/cli.rs`.
- Main app orchestration: `crates/gammaloop-api/src/lib.rs`
  (`OneShot::load`, `OneShot::run`).
- Python module entry: `crates/gammaloop-api/src/python.rs`
  (`#[pymodule(name = "_gammaloop")]`).
- Python package shim:
  `crates/gammaloop-api/python/gammaloop/__init__.py`.

=== 1.1 CLI architecture
<11-cli-architecture>
- The CLI is state-centric: every session operates on one active state
  folder and one active in-memory `RunHistory`.
- Startup is orchestrated by `OneShot`:
  - resolve the state folder from CLI / boot card / defaults
  - load or create the state
  - load persisted settings and run history
  - optionally apply a boot card before entering the REPL or running a
    one-shot subcommand
- Command execution is centralized in
  `crates/gammaloop-api/src/session.rs` (`CliSession`):
  - top-level commands
  - boot-card replay
  - nested `run` execution
  - command-history recording policy
  - command-block-definition mode
- The REPL layer in `crates/gammaloop-api/src/repl.rs` is state-aware:
  - prompt text reflects the active state settings
  - tab completion is driven from the active state, active command
    blocks, and clap argument metadata

=== 1.2 Python session architecture
<12-python-session-architecture>
- The `_gammaloop` module exposes a stateful `GammaLoopAPI`. One instance
  owns its `State`, `CLISettings`, `RunHistory`, default `RuntimeSettings`,
  and `CliSessionState`.
- Construction delegates startup to `StateLoadOption::load`, so Python and
  CLI sessions share state-folder, settings, boot-card, and logging setup.
- Command strings run through `CliSession`; structured sample evaluation
  uses the shared `EvaluateSamples` command. Python integration uses
  `api.run("integrate ...")`. The retained typed `GammaLoopAPI.integrate()`
  entrypoint remains disabled in the source. Its command-building helper covers
  slot selection, workspace/resume settings, and optional Gaussian-reference
  acceptance; signed and componentwise absolute integration-result classes
  remain registered, without a supported structured Python integration producer.
- `GammaLoopAPI` is a PyO3 `unsendable` class: each instance remains bound
  to its Python thread, while native work invoked by it may still use
  internal parallelism.
- `read_only_state` blocks writes inside the state directory. It does not
  make the session immutable: commands may still change in-memory state,
  settings, caches, observables, and history.

=== 2. Application State and Command Model
<2-application-state-and-command-model>
- Central mutable app state: `crates/gammaloop-api/src/state.rs`
  (`State`).
- Command history and run cards: `CommandHistory`, `RunHistory`.
- Command dispatch: `crates/gammaloop-api/src/commands/mod.rs` with
  concrete subcommands in `commands/*`.

The command model is stateful by design: commands mutate a long-lived
`State` that can be saved and resumed.

Selective loading filters processes and integrands before loading their heavy
payloads. The retained in-memory process IDs are dense, and loaded generation
reports are rebound to those IDs. This selection uses the same manifest and
version validation as a complete state load.

Generation preserves the original forward sides by default. The optional `--symmetrize-left-right-states` CP optimization remains a user assertion about the selected theory, process and coupling point, with warnings at generation and runtime warm-up. Models declare covariant cut multiplets, and generated integrands retain their physical event representatives. Regenerate saved processes and integrands after the phase/model changes; see #link("phase-conventions.typ#generation-options-and-generated-states")[the generation options and generated-state contract];.

=== 3. Domain Core (gammalooprs)
<3-domain-core-gammalooprs>
- Root module wiring: `crates/gammalooprs/src/lib.rs`.
- Initialization and shared symbol registries:
  `crates/gammalooprs/src/initialisation.rs`.
- Diagram generation and filtering: `crates/gammalooprs/src/feyngen`.
- Model and parameters: `crates/gammalooprs/src/model/mod.rs`.
- Graph domain: `crates/gammalooprs/src/graph/mod.rs` and submodules.
- Momentum routing and parameterization: `crates/gammalooprs/src/momentum`.
- CFF construction, numerator processing, and subtraction:
  `crates/gammalooprs/src/cff`, `crates/gammalooprs/src/numerator`, and
  `crates/gammalooprs/src/subtraction`.
- Process orchestration: `crates/gammalooprs/src/processes/mod.rs`,
  `crates/gammalooprs/src/processes/process.rs`.
- Amplitude and cross-section pipelines:
  - `crates/gammalooprs/src/processes/amplitude.rs`
  - `crates/gammalooprs/src/processes/cross_section.rs`
- Integrand abstraction and implementations:
  - `crates/gammalooprs/src/integrands/mod.rs`
  - `crates/gammalooprs/src/integrands/process/mod.rs`
  - `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`
  - `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`
- Integration engine: `crates/gammalooprs/src/integrate/mod.rs`.
- Global/runtime settings: `crates/gammalooprs/src/settings/mod.rs`,
  `crates/gammalooprs/src/settings/global.rs`,
  `crates/gammalooprs/src/settings/runtime.rs`.
- Differential event/observable pipeline:
  - `crates/gammalooprs/src/observables/events.rs`
  - `crates/gammalooprs/src/observables/mod.rs`
  - `crates/gammalooprs/src/observables/clustering/*`
  - `crates/gammalooprs/src/integrands/evaluation.rs`

=== 3.1 UV renormalization
<31-uv-renormalization>
UV generation defaults to the disconnected-capable hedge-poset backend
and also supports a legacy DAG-forest backend and a comparison mode.
Scheme policy remains on each Spinney while compute nodes store typed
local, integrated, and cut-dependent results. Four-dimensional
disconnected counterterms factorize over complete component
counterterms; three-dimensional disconnected terms replay
component-local operations from their common root so CFF structure is
not multiplied as though it were scalar. The maintained sign,
projection, marker, and backend-boundary invariants are documented in
#link("uv-renormalization.typ")[`uv-renormalization.typ`].

=== 3.2 CFF production and numerator-energy ownership
<cff-production-and-numerator-energy-ownership>
GammaLoop owns production graph/source construction, UV orchestration, exact source mapping, and evaluator preparation. The `three-dimensional-reps` crate owns the shared CFF algebra. The `3Drep` command and feature-gated eager evaluator are diagnostic tools, not production contracts: GammaLoop may prepare their inputs, factors, and expressions differently.

The shared `LinearEnergyExpr` stores exact `Rational` coefficients for indexed internal/external energies, the uniform scale and the constant term; `CFFVariant::prefactor` is also `Rational`. Arithmetic and cut handling retain that type until symbolic output converts it with `Atom::num`. Native rational serde/bincode support owns coefficient persistence; old Atom coefficient encodings are not a compatibility contract.

CFF capacities belong to independently sampled denominator occurrences. Physical sources use their EMR/source-edge identities; completed UV sources additionally accept typed canonical denominator classes, which are distinct from `EdgeIndex`. LMB coordinates certify routing and fixed affine carriers. They do not authorize redistributing a physical source's energy powers or combining contours.

Both direct local-3D modes first perform the complete loop-energy integration and build the complete/global CFF expression; the UV Taylor operators then act on that CFF expression. They use the same Taylor-transformed CFF bodies. Writing the generalized residue map as `{ k -> C_k }`, `explicit_orientation_sum_only=false` is `sum_k sigma(k) C_k`, with `sigma(k)` the one-hot selector for the complete residue-map key. The Taylor operator is applied independently to every keyed body and leaves that selector opaque. `explicit_orientation_sum_only=true` only replaces each selector by one and explicitly sums the same bodies. Neither direct mode reconstructs or projects completed local-4D Taylor structures.

`Integrands` retains flat, shared numerator definitions alongside its cut-keyed roots. A root calls `numerator_family(scope, coefficients...)`; its definition is an ordinary factorized tensor expression with scalar formal coefficients. For the current CFF maps these include the sign multiplying an on-shell energy and the integer coefficient multiplying the independent sampling scale `M`. The binding uses the complete `LinearEnergyExpr`, including rational loop-basis and external-energy coefficients and constant shifts. It therefore does not assume that a numerator depends only on energies inside its UV component, and does not implement an LTD consumer. Different complete maps remain distinct even when their production selector host is equal.

Root-only operations attach prefactors without multiplying the shared bodies. Semantic substitutions and derivatives visit both roots and definitions; exports and algebraic comparisons explicitly resolve the definitions. Saved amplitude data retains the definitions with its roots. Older saved-state layouts must be regenerated.

The projected local-4D route retains raw `Local4dCts`, recursive sectors and original provenance for subsequent outer Taylor operations. Its projection view normalizes completed hard roles zero and one into exact signed denominator classes. Each class includes the component domain, routing, mass and full polynomial; physical/soft provenance, frozen localizers and cograph bindings remain explicit. The class algebra contains powered denominators and a factorized numerator, with no CFF capacities or sampling conventions. This is the boundary a future LTD consumer can use directly.

Local Taylor construction retains spinor products, chains and traces after metric simplification. Numerical tensor execution contracts them after residue mapping; the integrated-CT preparation separately performs the required analytic Dirac algebra on its copy before Vakint. An independent color pass reduces color contractions before shorthand traces are opened, preventing equal color spectators from acquiring distinct dummy-index copies while keeping kinematic coefficients factorized. This avoids expanding a local trace into scalar contractions before canonical merging and energy sampling. The local 3D, infrared and local 4D Taylor calls use the existing `NumeratorAtomExt` to retain an outer product's factors that do not depend on the expansion variable. The native series engine expands the complete dependent product, including Laurent poles, and the independent product multiplies its truncated Atom once. Outer sums retain native coefficient collection, preserving literal zeros exposed by that collection. This boundary does not move a dependent mapped numerator outside Taylor or collect distinct denominators together. For a retained numerator family, direct-3D Taylor expansion computes coefficient bodies once with the energy-map coefficients left formal. Finite absolute-order Laurent probes give lower bounds for both bodies and complete products, including denominator poles. The resulting coefficient families are specialized only by their residue arguments. Independent tensor factors stay outside coefficient collection. Known multilinear dot and projector arguments are normalized before the native series call; temporary linear heads never escape that operation. Specialization must be polynomial in the formal map coefficients, so it can cancel leading terms without introducing unaccounted Laurent poles. Analytic spin expansion temporarily aliases completed, spin-independent scalar products, certified by their expanded tensor structure. Their internally bound indices must not be reused when a powered product is materialized. Open compact tensor contractions remain visible to full shorthand expansion; the existing alias owner restores scalar products afterwards. Tensor execution evaluates odd powers as paired contractions times the remaining base, preserving its free indices and the requested exponent. Both local routes simplify the reduced numerator's color algebra before Taylor construction and the cograph's color algebra before residue mapping. The direct 3D reduced-kernel call and both final-cograph calls disable `ColorSimplifySettings::simplify_non_color`, retaining the factorized momentum coefficient while reducing the collected color payload. Open color indices remain explicit; the final color pass contracts indices closed by attached UV terms and projectors. When fundamental-dimension invariant substitution is requested, the color simplifier also applies it after each local rewrite and before tensor collection. Resolved scalar Casimirs and indices therefore do not cause collection over a factorized residue sum; newly produced invariants follow the same rule. Raw graph storage and parse-time validation are unchanged.

Tensor collection temporarily aliases maximal unselected composite coefficients before polynomial grouping. Selected tensors retain their complete slots and payloads; structurally equal coefficients share an alias, with collision checks against input symbols. The existing Symbolica alias owner restores definitions before tensor-group callbacks run. This preserves factorized powers and sums and avoids statistical zero tests on large momentum coefficients. Collection does not promise polynomial simplification of those opaque coefficients; structural cancellations and tensor-algebra zeros still apply. Aliases live only for the collection call and are never serialized or retained in generation caches. Chain composition uses this same collector before applying its existing composition and normalization rules, including when color simplification encounters chain factors in a mapped numerator. Positive integer powers of collected tensors stay compressed inside the complete payload passed to callbacks. The color simplifier runs its fixed-point rewrite loop on these collected payloads, keeping unrelated momentum coefficients outside repeated chain traversal. Untyped color heads follow the same selection. With the default settings, generic non-color chain/trace and metric rewrites reach their own fixed point before and after color collection, without making non-color trace sums into polynomial variables. Calls that disable non-color simplification defer these rewrites to the later normalization boundary. Chain-projector structure inference combines the exposed slots of its factor sequence when no direct slot or index-bundle arguments are present. Enabled chain and trace materialization expands the projector through the existing projector owner before constructing links, so symmetric, antisymmetric and cyclic sequences retain their external indices.

Final residue assembly validates the complete cut-key shape of every summand, then uses Symbolica's native bulk addition once per cut order. Projected sectors and direct-route selector materialization retain their already mapped branch values until this merge. This avoids repeatedly copying the growing numerator; it merges existing top-level sums without distributing product factors or constructing common denominators. After all selected forests are assembled, final factor collection protects complete functions and powers. Inverse-containing factors and selectors stay local to each occurrence; tensor sums under products retain their original contraction boundary. The collector therefore recovers common numerator factors without clearing denominators, sharing branch guards, or pulling a summed vector across independent contractions. Repeated collection stops at a fixed point.

Vakint's pole-order query uses absolute series depth zero, which includes every pole without searching indefinitely for a nonzero leading term in an exactly cancelling coefficient. The complete coefficient remains available to the separate finite/evanescent Laurent expansion; that expansion's depth is unchanged.

Tensor-network edge joins use Linnet's existing dense edge swap, which updates the moved edges through their owning half-edges. Each join therefore avoids a scan over the growing graph while preserving its edge payloads and merge flow. Operator merging uses its existing disjoint union groups to mark newly internal half-edges directly from each group's node crowns. One shared deletion set preserves original self-loops, dangling slots and edges between groups without allocating and combining graph-sized masks for each operator island.

Ready-operation batches likewise retain sorted visible half-edge lists per operation. Disjointness checks and the optional batch union visit these lists; counts use their lengths. Node identification consumes the same ordered list and its exact membership predicate, preserving the dense subset's representative node and newly internal edge marks, including partially hidden edges. Sequential, parallel and partial-rewrite strategies use this common payload. Only extraction materializes a dense subset for an individual operation, with the complete graph extent. Retained batch membership therefore scales with touched incidence rather than the product of graph size and ready-operation count. Root alignment uses the existing traversal's discovery order to collect one root per visited node, without converting that traversal to child vectors. Requests for disjoint node crowns commute. The traversal is released before the graph's node-store conversion, which still normalizes sibling ordering.

Rational-shell extraction groups equal denominator multisets through sums, products and powers while leaving denominator-free numerator subtrees opaque. It adds powers and numerators, prunes exact zeros, and never constructs a global common denominator. Positive typed denominator factors remain indivisible numerator blocks with their complete polynomial and signed binding. Compatible sectors merge before reconstruction. Exact cancellation at this boundary returns the existing typed zero with all cut-order outputs preserved, before any energy map or orientation is required. Frozen domains and localizers remain part of the compatibility test; opposite terms in different domains do not cancel. After each independent component contour, states with equal remaining requests and bindings share numerators when their carriers agree, or share carriers when their numerators agree. Grouping repeats until the state count stops decreasing, because one sum can create an equality for the next grouping pass.

Every algebra bucket retains one deterministic physical source witness. The existing source builder inherits original incidence, contracts omitted edges, and adds serial copies for raised occurrences. Only certified serial paths or pure cycles permit incidence contraction of equal channels. Nonadjacent equal channels can merge algebraically while retaining distinct edges in that witness. Exact signatures validate routing and its sign; they never infer endpoints. The source separately checks rank, measure, boundary attachments and contour. The class numerator is certified in one neutral momentum frame, independently of denominator routing, mass, multiplicity and domain. The signed occurrence lift obeys `H=hR`, `P=rR`, hence `H^0=hrP^0`.

Completed hard numerator factors may use any certified occurrence of their class, including occurrences inherited from different physical owners. Equal channels with unequal masses or prescriptions remain distinct pools. A hard carrier without a surviving class uses an exact fixed affine lift through retained base occurrences and external shifts; it creates no new pole family. Physical source ownership and soft crown provenance retain their restrictions. Non-vacuum exact sources retain pure-external boundaries as explicit source-crown hedges. Future on-shell insertions such as `(m,0,0,0)` require an explicit fixed boundary payload rather than inferred attachments.

The immutable CFF assignment is constructed compositionally. Sums take the componentwise maximum, products and registered multilinear arguments add degrees, and positive integer powers retain compressed repetitions. A deterministic checkpoint with a doubling interval detects assignment cycles after a transient while retaining one offset context, including packed tails beyond 64 occurrences. The inert cyclic, symmetric and antisymmetric tensor projectors use the same slotwise degree rule while keeping their symbolic heads intact. Freely assignable leaves use cyclic offsets within each class; opaque blocks use deterministic least-loaded placement. An ordinary affine request admits one candidate. A nonlinear request admits at most three distinct valid candidates: baseline, packed unavoidable excess, then reversal of the best-scored placement (rotation only replaces a duplicate). Each challenger changes one eligible occurrence pool, with no Cartesian search across pools.

Selection uses native `generated.expression.orientations.len()` before surface conversion and host selection. Maximum rank, descending rank envelope and stable candidate order break ties. The cache key preserves exact ordered capacities; identical physical energies do not authorize permuting loads. This bounded selection is not an exhaustive minimum. The sunset regression compares complete signed contours against independent residues: at energies `(2,3,7)`, the degree-five and degree-seven values are `1213/16387080192` and `-365/47775744`. It covers all 126 degree-five distributions over the five equivalent occurrences and the allocator's degree-five/seven proposals. It imposes no native row count or proposal order. All losing generation belongs to selection cost. Soft routing retains its existing exact-basis/Pareto frontier for off-shell proposals, including every fixed external shift. Its at-most-three native generation trials do not bound that separate frontier's preparation cost. Root-expression reuse retains its established capacities.

Generation profiling includes raw physical-degree reports and the complete outer routing selection interval. Raw CFF requests separate source reconstruction, degree/preparation, native generation and postprocessing, and selection records the winner explicitly. This keeps winning generation separate while charging every losing proposal, its preparation and destruction to dispatch. Full local UV accounting uses completed nonroot forest nodes; child projection timers alone do not include the outer routing and assembly boundary.

One nonserialized `Local4dProjectionContext` spans each graph's complete UV computation in both forest orchestrators. Its deterministic LRU caches retain canonical sectors, owned source analysis and immutable numerator templates together (48 MiB, 4,096 entries), winning CFF payloads and contender counts (64 MiB, 4,096 keys), and component-local mapped subtrees (16 MiB, 16,384 entries). Owned keys, containers and symbolic payloads are charged; oversized entries bypass retention. These are retained-payload limits, not RSS limits. Clearing or evicting entries changes computation cost only. Source preparation owns its parsed incidence and signed namespace maps, so a warm hit skips reconstruction without retaining a borrowed graph adapter. The context belongs to one graph computation and cannot be reused for another graph; its keys describe that graph's requests rather than a global graph registry. Projection reuse tests compare complete coefficients with cold, warm, disabled and evicted retention, while checking only the configured byte ceiling. Hit counts, prepared entry counts and traversal order remain diagnostic observations.

`PlannedExactSourceNumerator` owns the assignment, prepared expression and exact signed mapping context behind an immutable `Arc`. A cache miss prepares the same generic template. Capture-free loop/occurrence parameters remain symbolic until a residue row arrives; each subtree records its exact dependencies. Mapping converts required energies once per row and memoizes only the relevant sample tuple. Constants are reused unchanged. Zero contact samples and the order of inactive-energy elimination are preserved. Mutable row caches live outside the immutable template and are cleared after each component wave. No graph-specific tables, global caches or serialized caches are involved.

Projected local-4D assembly maps one ordinary numerator template per compatible component family. Subsequent components see that template with earlier residue coefficients still formal; joint rows append only scalar arguments and CFF carriers. Families merge only when their remaining denominator powers, ordinary body and ordered formal signature agree. The final outer-CFF preparation keeps the ordinary numerator and scalar carrier visible as separate factors of the same certified routing product, then binds both with the same complete outer map. Derivatives can couple components, so this representation generally retains a joint residue sum. It does not replace that sum by a product of independent sums without a dependency proof, or copy the complete numerator into every row.

Each candidate prepares or reuses its exact immutable template before native CFF generation. Preparation certifies the occurrence diagonal in a formal spatial/temporal frame, including fixed external four-vectors. The selected payload carries this same certified template into row mapping. Profiling separates winning template construction from discarded candidate construction; all discarded preparation and cache work count toward dispatch overhead. Its only distribution is the finite projector action on an affine energy sum; numerator products and powers remain intact. A positive denominator block with a pinched hard carrier keeps its whole polynomial and an exact source-owned affine-rewrite witness. Only that witness permits joint fixed dependencies in the existing planned factor; uncertified mixed ownership and soft provenance keep their original rules.

The full reconstruction invariant, sign argument for `D(Q)=D(-Q)`, and worked fixtures are documented in #link("exact-powered-denominator-cff-lifting.typ")[`exact-powered-denominator-cff-lifting.typ`];.

The production numerator remains factorized. Degree analysis traverses its factors without expanding them, each UV step attaches only newly owned factors, and final assembly attaches outside and global factors exactly once. When one outer CFF serves independently evaluated residue branches, its capacity is the per-edge maximum of their separately analyzed factorized products, including the common remaining numerator. Equal branch bodies need only one rank analysis; opposite bodies remain independent and cannot cancel each other's capacity. The analysis streams degree maps without constructing a tagged symbolic sum. Initial-cut and tree-edge energies remain outside this internal CFF capacity domain; stored production-root CFF reuse retains its existing boundary. For higher power projection, the term parser splits only the completed expression's outer Taylor sum when its addends carry separate denominator topologies. Nested numerator sums and positive typed `GS.den` factors remain factorized and are never cancelled against denominator occurrences upstream; generalized CFF owns the resulting pinches and lower sectors. For higher powers, interpolation may replace an EMR energy by `a*M`, where `a` is a signed integer and `M` is the common auxiliary CFF numerator-sampling scale. This is an EMR substitution, never an LMB rewrite. Production evaluator parameter lists always include `M`, and runtime settings require `M != 0` for every process. The physical result is invariant under changing its nonzero value.

The signed sampling coefficient is neither a physical edge-direction sign nor the runtime residue-map-key selector `sigma(map_id)`. Pole signs enter the generation-time affine map, after which an orientation entry stores `a` in `LinearEnergyExpr::uniform_scale_coeff` and its interpolation weight stores the compensating inverse power in `CFFVariant::uniform_scale_power`. Multiple numerator maps may therefore share the same coarse orientation. Their complete loop/edge energy maps and distinct `numerator_map_index` values remain attached to the factorized numerator. Each such entry has its own complete map key and therefore its own `sigma(map_id)`; entries must never be merged merely because their physical edge directions agree.

Integrated finite UV terms retain their exact source-local EMR maps. Production map-key selectors partition complete generalized residue-map entries but do not replace those maps. In orientation-local direct-3D generation, shrinking a UV subgraph preserves the orientation selected on every surviving outer edge. Among the cut-valid full orientations which extend that outer assignment, one deterministic inner orientation is selected to host the complete integrated finite counterterm; the other inner extensions receive none, so summing orientations counts the addback exactly once. In explicit-sum direct-3D generation the same reduced residue is kept once without a selector. Projected local-4D counterterms are likewise selector-free: each completed exact source owns its full source-local orientation sum, independently of production-orientation IDs. After that sum, a projected sector stores only its coefficient and frozen localizing factor. Its child map lengths need not agree across Taylor topologies. Actual outer branches keep their own maps until the coefficient and remaining numerator have been mapped, then sum by cut order. Typed zeros retain the allowed cut keys without requiring a map; direct and integrated localization retain their own production-host rules. The empty UV forest is the ordinary factorized production root in both local-UV routes; the expanded-4D setting changes only proper, nonempty UV nodes.

Final assembly uses one `Integrands` map of cut-indexed factorized expressions. UV markers and tensor replacements act on those expressions before the ordinary evaluator preprocessing; there is no parallel deferred-body state. Direct-3D sequential forest construction composes the same local and integrated replay operations used by disconnected forests. Integrated localization registers its source surfaces before either Taylor branch reads the graph, and its normalized localizing factor remains inert under subsequent Taylor operations. Evaluator preparation skips residue-selector parametrization when the selector symbol is absent, and skips orientation collection when neither theta nor IF occurs. Expressions with conditionals retain the existing guarded-branch rewrites so inactive residue-local inverses are never evaluated.

The shared CFF core also returns its connected-loop and pure duplicate-denominator global sign as typed metadata. GammaLoop consumes that bridge exactly once for root, reduced, and exact production CFF sources, cancelling the shared-core-local uniform convention before the physical Minkowski measure `i^L/(2*pi)^(3L)` is applied. Integrated UV addbacks use the same convention, with additional Vakint normalization `1`. The NLO acceptance layer independently generates orientation-local direct 3D, explicit-sum direct 3D, and projected local 4D with local and integrated UV and threshold counterterms. It compares complete GL0/GL2 values at a common native- Arb point in all three routes; the fast Monte Carlo tests integrate explicit-sum 3D. DD acceptance checks the inclusive `(alpha_s/pi) * LO` correction, graphwise UV-mass and localization-scale independence, cancellation of total renormalization-scale dependence, opposite GL0/GL2 squared-scale logarithms, and the physical and projected EMR energy bounds. TT acceptance uses the fully-MSbar scheme, without on-shell counterterms.

Direct-photon benchmarks use the off-shell spin projector `-g^(mu nu)` and no picobarn conversion. The inclusive lepton-process targets instead use the Eq. (7.1) normalization `2(4 pi alpha)/(3 Ecm^3)` and the conversion to picobarns; individual lepton-process graph components are not assigned the unconverted published photon targets. Acceptances require positive real LO and the signed real NLO correction, checking the imaginary component separately. Diagram contributions can have either real sign. Current validation results and measured timings belong in the accompanying test evidence and PR.

The raw forward graph already contains inverse-process UFO vertices with their Hermitian-partner spin structures and couplings. An additional RHS numerator adjoint would conjugate those structures a second time. Marking vertices and virtual propagators, together with the reversed RHS virtual contours, instead gives `(-1)^C_R`, with boundary hairs retained in the RHS component count. The complete LU residue factor is `2*pi*i*(-1)^C_R`, hence `-2*pi*i` for connected RHS, while retaining the cut propagator numerators exactly once. Bare, local and integrated UV, and threshold branches use this same cut-group factor. Runtime subtracts threshold helpers with left `-i*pi` and right `+i*pi` integrated coefficients; iterated terms use their product. See #link("phase-conventions.typ")[phase conventions and the independent conjugation audit] for the cut-line `i` cancellation, intrinsic complex CKM phases, and the unsupported complex-mass-scheme boundary. There are no conjugation modes.

The scalar local-equivalence matrix is generated from the scalar model rather than from hand-built graph data. Its lanes without an additional numerator retain the full UFO Feynman rules; companion probes use only Feynman-rule-local edge factors, and there is no graph-specific production branch. The matrix enables local UV, integrated UV and threshold counterterms while comparing all three local-UV routes, including native-Arb checks. `just test_LU_scalar_xs` includes the slow cases, uses release compilation by default and stops on the first failure. Nonzero route comparisons require finite values and precision-scaled relative agreement. Only an independent exact source-zero certificate permits the separate absolute bound; small magnitude alone does not qualify. The curated suite includes the six base scalar graphs GL00, GL02, GL04, GL08, GL09 and GL24. Each profiles direct local3D separately for every complete residue-map key and projected local4D after the complete residue sum.

UV profiling defaults to `only-divergent`: every expected cycle union with DOD \>= 0 is tested using the generation LMB when suitable, otherwise the first suitable basis in the deterministically sorted complete LMB list. The exhaustive `all` mode is opt-in. Amplitude and LU inputs share this behavior, with graph and Cutkosky-cut selectors for LU profiling and a colored final failure summary. Per-key profiling is defined for orientation-parametric, localized direct local3D. Selector-free explicit-sum direct local3D and projected local4D are summed representations and reject that request.

=== 3.3 Tensor-network contraction order
<tensor-network-contraction-order>
Generic metric normalization retains its ordered vector, redundant-metric, vector-power, metric-power and trace passes. The first two passes traverse the factorized expression from the root downward and use conservative symbol-tag checks before applying their unchanged rules at the current node. Ordinary momentum components and scalar dots avoid unnecessary wildcard matching. A matched replacement keeps the original descendant-skipping behavior, and the redundant-metric pass still runs to a fixed point. In particular, metric-power normalization continues to precede tracing.

Evaluator construction normalizes its factorized input once, then parses each existing top-level summand into an independent Spenso network. This bounds network preparation to the current summand. Products, powers and nested sums retain their grouping. Parsing and contraction use `ParamTensor<ShadowedStructure<Aind>>` throughout, with symbolic `Atom` components. The HEP tensor library supplies `ParamTensor<ExplicitKey<Aind>>` entries, including exact metric and identity components. Numerical evaluation remains a later evaluator step.

Shared numerator families are lowered before this outer contraction. Preparation first groups complete contractions of family calls with their surrounding tensor factors, including trace and chain contexts. Closed additive scopes are prepared from the inside out and remain factors in surrounding products. The existing tensor collector exposes residue sums over opaque family calls while preserving inline numerator sums, scalar coefficients, independent closed products and conditional guards. Only a closed tensor interface becomes a shared joint body. Each call occurrence receives fresh scalar formals and private contraction indices; exposed ports retain their identity. Existing chains and traces own their implicit endpoint markers: chainification keeps those scopes opaque, and joint preparation rejects substitutions that would introduce nested binders. A structural template reuses this body across physical argument rows. Parameter-dependent tensor ports and formals captured by existing aliases retain their original family boundary. Only referenced bodies are then contracted, once per definition in an evaluator stack, and any remaining open tensor indices are represented by temporary tensor calls. The outer contraction closes those indices and selects the corresponding scalar component functions. Only generated functions reachable from the final roots and scalar aliases are registered; outer projectors can discard otherwise populated components. Component and scalar-alias signatures retain only their required formal arguments, including transitive alias dependencies; physical family signatures still describe the complete energy map. Known sparse zero components vanish, while unknown families and invalid index dimensions remain errors. All evaluator modes receive the resulting function-map definitions. The pinned Symbolica translator inlines these scalar functions when generating instructions, so a compact input does not by itself guarantee fewer runtime operations. Closed shared bodies avoid enumerating their open tensor components; generation and runtime effects still require measurement.

Before parsing a summand's complete contraction, evaluator preparation may materialize an existing additive tensor factor of its root product. It inspects the parsed factor's actual scalar and tensor entries and admits only Gaussian integers (integer real and imaginary parts), with a nonempty self-dual exposed boundary. Exact Gaussian-integer components retain symbolic arithmetic at arbitrary integer size. Floating coefficients decline this preparation and retain their original contraction order. Conservative Cartesian-capacity and component-size estimates reuse the 16 MiB eager tensor-sum budget. Explicit lazy tensor sums skip this preparation. Qualifying factors use the same execution and contraction strategies as the complete network, retaining their symbolic components and index ownership. Existing network multiplication joins them to the grouped residual before scalar aliasing and sum-boundary closure. When no factor qualifies, preparation parses the original summand directly. Taylor and other raw symbolic callers do not use this finite component preparation.

Each network aliases large scalar references and contracts its tensor products, then carries the root and definitions in Symbolica's existing `AliasedAtom`. Selectors stay visible outside these scalar definitions. Distinct network, summand, output and concrete-orientation scopes prevent handle collisions; fresh symbols also avoid registered alternate names. Handle lookup avoids hashing subtrees larger than any key; renaming still visits their children. Explicit alias resolution retains its exact fixed point for nested or newly exposed handles, avoiding unrelated hashes and copies of alias-free results. Within a network, bulk Atom addition accepts every leaf whose existing scalar conversion succeeds, including rank-zero tensors and lazy sums with scales. The conversion preserves alias handles and rejects open tensors before the store changes. Owned scalar results move into the existing bulk or streaming sum, avoiding repeated pairwise addition of large closed-tensor results. Scalar roots and their definitions are combined before global evaluator optimization; open tensors, including open zero tensors, remain invalid scalar outputs. Optional evaluator variants borrow these scalars before the final parametric evaluator takes ownership. Root substitutions also apply to alias bodies. Concrete residue-map selection precedes physical orientation signs and prunes unused definitions, keeping inactive singular branches out of evaluation. After tensor contraction, scalar sums with the same complete residue-map condition share one lazy guard. Their bodies remain inside that guard so numeric construction can share calculations without evaluating an inactive contribution. The actual Symbolica evaluator registers the retained definitions in its function map. When source storage is enabled, stored roots and complete function-map entries together form the source representation, including for archive reloads. Selector-free inputs and empty function-map replacement lists retain that ownership without copying the expression. Once symbolic programs are built, unstored source expressions are released before dual and numeric programs are constructed. Profiling separates expression preparation, Symbolica construction and numeric-program conversion. At this finite component boundary, a ready tensor whose exposed indices all contract with a pending tensor sum is attached to each immediate sum branch first. One expression traversal selects all eligible disjoint outermost sums, including siblings within a product. Small product wrappers copy only the closing leaf references and their exact slot order, flow and internal traces. Native arms and sum shells stay in place; stable half-edge identifiers bind the consumed inputs while residual seams remain untouched. The wrappers are appended together, followed by one deletion and operator merge per wave. This reduces the sum's open tensor rank before large component expressions are constructed, without repeated extraction of a whole sum or its arms. The tensor store is shared by reference and scalar spectators stay outside; sums are not multiplied through other sums, powers or opaque functions. Repeating the transformation decreases an existing sum's exposed rank. Eligibility requires self-dual exposed sum slots, so shared edge descriptors retain their exact endpoint meaning; internal dual contractions and spectator incidence are preserved. It changes the contraction graph, leaving the input Atom, raw Taylor algebra and selected energy assignments untouched. Symbolic networks used to reconstruct numerator Atoms do not invoke this preparation.

Completed contraction waves reclaim tensors no longer referenced by the whole network graph. The owned store moves surviving entries and remaps every tensor leaf variant, including shared references and scaled sums. Tensor indices are local execution handles and may change between waves. Reclamation occurs only after replacement leaves are installed and deferred node identifications are finished; parallel workers first return all additions from their borrowed overlays. Extracted subgraphs cannot determine whole-store liveness. Scalar indices remain stable, preserving alias handles embedded inside expressions and the original definitions used during final restoration. No symbolic expression is rewritten by this storage compaction.

The default `intermediate_cost` preset selects `MinIntermediateCost`, a configuration of the existing `ContractionStrategy` implementation. It changes pair selection, preserving the scalar, trace, library, disconnected-product and final cleanup rules.

The greedy score first estimates copied symbolic payload, then symbolic normalization work, before using sparse support, tensor volume and rank as later criteria. For a candidate pair, let `J` count the estimated nonzero entry products with matching contracted coordinates. Let `b` and `t` be each operand's mean entry bytes and top-level terms, rounded up. The two leading scores are `J * (b_left + b_right)` and `J * (b_left * t_right + b_right * t_left)`. Sparse support supplies the join count without multiplying coefficients. Otherwise the fallback uses the Cartesian entry count, capped by output and contracted-coordinate volume for single tensors. That cap is not applied to lazy sums, whose terms may share coordinates. The existing bounded output-support join and deterministic operand-order ties remain in use. The bounded join counts output coordinates through borrowed slices of its immutable support groups, avoiding per-product coordinate copies while merging equal output coordinates from distinct contracted groups by their contents.

These scores estimate work before cancellations. They are not physical memory bounds or a prediction of the globally best contraction sequence. They use existing structure and scalar-size profiles, with saturating arithmetic; the planner performs no speculative contractions or numerator expansion. A higher rank intermediate can therefore be preferred when its estimated cost is lower. Performance claims require identical-input measurements of both time and memory.

The previous sparse-aware preset and all other presets remain explicit choices. To restore the preceding production configuration in a run card, use:

```toml
[cli_settings.global.generation.evaluator]
spenso_execution_mode = ["Sequential", "MinResultRank"]
tensor_network_contraction_order = "sparse_atom_aware"
```

Within a running session, switch the preset with the existing command `set global kv global.generation.evaluator.tensor_network_contraction_order=sparse_atom_aware`. Use `intermediate_cost` to select the default again. The preset applies to `ContractionMode::MinResultRank` with the supported sequential executors; explicit `SmallestDegree` selection is unchanged. The existing `bnl_evaluator_atom_mwe` diagnostic exposes the same choices through `--contraction-order`, using hyphenated names such as `intermediate-cost`.

== Lifecycle and Data Flow
<lifecycle-and-data-flow>
=== 1. Startup
<1-startup>
+ CLI parses command input (`OneShot::parse_env_with_capture`).
+ Initialization runs (`initialise()`), installs hooks and symbol
  registries.
+ State is loaded from state folder if available (`State::load`),
  otherwise default state is created.
+ CLI settings/runtime defaults are loaded and synchronized with tracing
  filters.

=== 2. Process Generation Flow
<2-process-generation-flow>
+ `generate` command builds `ProcessDefinition` (from syntax or graph
  import).
+ `State::generate_integrand(s)` creates a generation thread pool.
+ For generated processes, `feyngen::DiagramGenerator` constructs graphs
  from model vertex rules, filters them, and performs topology- and optional
  numerator-aware grouping. Graph import supplies that boundary directly.
+ `ProcessList::preprocess` delegates to the amplitude or cross-section
  pipeline. Those graph-level stages generate CFF/cut surfaces,
  loop-momentum bases and parametric integrand data, plus the configured
  threshold- and UV-subtraction data.
+ `ProcessList::generate_integrands` packages the resulting graph
  collections, runtime settings, and evaluators into `ProcessIntegrand`
  instances.
+ Optional compile/export steps persist compiled evaluator artifacts and
  DOT/standalone outputs.
+ Each generated integrand now embeds its frozen f64 backend choice in
  `integrand.bin`:
  - `eager`
  - `symjit`
  - external `c++` / `assembly` plus the external compile options used
+ Runtime-only evaluator backends are then activated from that frozen
  metadata:
  - eager uses the saved eager evaluator directly
  - symjit is rebuilt after generation/load from the saved Symbolica
    evaluator
  - complete external compiled artifacts are loaded while the saved
    state is activated
  - missing external artifacts leave the frozen backend metadata
    unchanged but activate the portable eager evaluator for the current
    session
  - if external loading fails and startup globals explicitly opt into
    symjit, GammaLoop falls back to symjit for that integrand and logs
    it

SymJIT rebuilds the saved Symbolica evaluator with its extra common-expression cache disabled: SymJIT 2.21 can reuse function results across inactive branches or discard stores still needed by later calls. Symbolica's numeric program keeps its existing sharing; requested optimization levels remain capped at O2.

=== 3. Evaluation and Integration Flow
<3-evaluation-and-integration-flow>
+ Commands (`inspect`, `evaluate`, `integrate`) resolve process +
  integrand references.
+ Integrand is warmed up (`ProcessIntegrand::warm_up`) to initialize
  rotations and caches.
+ Sampling path parameterizes points and evaluates graph terms.
+ Stability checks may escalate precision (`f64 -> f128 -> arbitrary`)
  and rotate kinematics.
+ Process graph evaluation returns a rich `GraphEvaluationResult<T>`
  rather than only a complex weight. This carries:
  - the graph contribution
  - grouped generated events
  - event-processing timing
  - generated / accepted event counts
+ A shared process-layer prepared-event helper is used by both
  cross-section and amplitude graph evaluation:
  - it decides whether an event must be built for the current rotation
  - runs selectors immediately on that candidate event
  - retains only accepted identity-rotation events when buffering is
    needed
  - reports generated / accepted event counts and event-processing timing
  - selector failures zero the local contribution, so rejected events do
    not survive into the final retained result
+ Stability selection checks the signed graph weight, retained absolute
  contributions, and Gaussian reference moments independently when present.
  All required checks must pass; the retained branch also carries the final
  grouped event payload.
+ The final `EvaluationResult` contains:
  - the stable `integrand_result`, including each map's Jacobian and channel
    partition factor, combined with the physical value at native precision
  - a unit unapplied `parameterization_jacobian` for x-space evaluation
    (`None` for direct momentum-space evaluation); the actual map Jacobian
    remains on the bridge's map evaluation for diagnostics (the selected
    momentum sample carries the combined map/partition factor)
  - the separate `integrator_weight`, i.e. the Monte Carlo/grid weight
    only
  - grouped accepted events
  - event-processing timing in `evaluation_metadata`
  - generated / accepted event counts in `evaluation_metadata`
  - ordered per-level stability results, each with relative-accuracy and
    total time spent in that stability level
  - evaluation metadata
  Native reconstruction and stability selection return the same retained
  precision-tagged result to both API routes. The ordinary route checks final
  f64 representability of contributions and event weights; the precise Rust
  route keeps the native result without a preliminary f64 reporting pass.
+ Observable filling happens only from the final stable
  `EvaluationResult` payload. Unstable branches do not contribute events
  or observables.
+ `havana_integrate` runs iterative Monte Carlo updates, merges
  worker-local observable accumulators, and writes integration
  artifacts.
+ Interactive integration status uses a renderer-neutral snapshot and a
  separately owned terminal event loop, as documented in the maintained
  #link("ratatui-integration-dashboard.typ")[integration dashboard architecture].

Every evaluator receives the auxiliary numerator sampling scale `M` as an input, defaulting to one. Warm-up rejects a zero runtime scale for every amplitude or cross-section integrand, including evaluators that do not use `M`. No source-expression scan or serialized usage flag is needed. This enforces the EMR-only `a*M` contract above; `M` is not an LMB coordinate.

Binary64 reporting may round exponentially suppressed, fully weighted components
to zero. Numerical factors are combined before this rounding so a compensating
factor cannot hide a meaningful contribution. For an averaged stability result,
the mean absolute fully weighted probe value can establish underflow below the
smallest normal binary64 value. This bound stays separate for each component and
observable; an imaginary scale cannot excuse a meaningful real discrepancy.
Overflow and invalid operations remain failures.

With `stability.check_on_norm`, amplitudes compare full complex magnitudes even
when integration selects only Re or Im: a fixed-helicity amplitude can acquire
a global phase under spatial rotation. A private sum of channel-point complex
norms separately detects instability hidden by signed channel cancellation.
Common-point physical terms are summed before taking each norm. The public
componentwise absolute result remains the primary probe's value. Cross sections
retain checks of the selected integration components, and Gaussian-reference
values and moments retain independent checks.

=== 3.1 Differential event-processing runtime
<31-differential-event-processing-runtime>
The differential pipeline is split into three layers:

- Declarative settings in `RuntimeSettings`:
  - `quantities`
  - `selectors`
  - `observables`
- A compiled runtime-only cache:
  - `EventProcessingRuntime`
  - owned by the process integrand
  - built in `warm_up()`
  - invalidated when runtime settings are mutated
- Per-sample / per-integration results:
  - `GraphEvaluationResult<T>`
  - `EvaluationResult`
  - observable accumulator state / snapshots

The runtime-only cache is intentionally not serialized in binary state
dumps. It is wrapped in a no-op `Encode` / default-empty `Decode`
container and is always rebuilt from settings at warm-up time.

The same pattern is now also used for evaluator execution backends:

- per-integrand frozen backend metadata lives inside `integrand.bin`
- loaded external `.so` evaluators are runtime-only
- symjit evaluators are runtime-only
- the saved portable representation remains centered on eager Symbolica
  evaluators

=== Canonical sampling channels
<canonical-sampling-channels>
`SamplingChannelCatalogue` supplies the one channel-ID domain used by grids, summed evaluation, explicit channel evaluation and event metadata. Generated LMB basis IDs describe graph routing only. Both amplitudes and cross sections use `SamplingChannelBridge` to map each selected unit-cube point into the parent frame and evaluate the partition at that raw point. Summed and Monte-Carlo evaluation apply `J_c w_c` once per graph result and event at native precision, before reporting, and retain a separate grid probability. Stability rotations act on the resulting mapped point and external frame together.

Both execution modes prepare each channel on the graph-group master and reuse that point for every group member. The absolute monitor first sums all physical cuts, counterterms, orientations and group members at that point, then takes the componentwise absolute value. Explicitly summed channels add these positive contributions before one outer-cube statistics update, retaining channel covariance. The signed channel sum is not passed through `abs` afterward. Graph groups remain separate integration domains; if orientations themselves are sampled, the absolute monitor describes that sampled-orientation domain, not the absolute value of an orientation sum. GL638 studies sum orientations. Native, precise and Python evaluation outputs expose the separate `absolute_integrand_result` payload before the outer grid weight. Integration workspace version 3 rejects older checkpoints whose absolute moments used the previous convention, requiring a fresh run.

The existing graph sampling setup retains a nonserialized canonical catalogue, compiled programs and native Double/Quad/Arb bridges. Process warmup constructs the configured precisions transactionally after numeric masses/externals are ready; mutable settings invalidate it and evaluation then requires warmup. Fresh bridge construction remains separate for inspection with explicit inputs. Workers clone compiled eager score buffers into separate mutexes, while each point borrows its maps and scores without cloning programs. External improvement clears old numeric caches before rebuilding from directly edited inputs.

Standalone `phase_space(cut(...))` maps use the real graph\'s cut equation, warmup masses and fixed external data through the shared implicit radial kernel. Their fixed directions use the existing LU ray: route signed velocities before radial scaling so large host-null master components cancel before multiplication. They retain the auxiliary raw radial variable. Conditional cut/left/right maps now share the Symbolica-resolved ordered block plan and one geometry registry. The registry key includes qualified target, native parent, active edges and ordered prerequisites. Exact routing rejects any unsampled dependency or side displacement that would change the host cut. The existing embedding prepares the host\'s native LU root from actual preceding coordinates, then maps physical active coordinates back to raw coordinates with the complete affine shift and determinant. Every foreign inverse prepares its own supplied raw point. The complete native map is routed to the master frame once. The combined generated fixture gates pass; general nested mixed compositions remain rejected. Explicit pairs of original graph surfaces can share one active three-dimensional joint block under a physical cut. Binding and adoption resolve the same original equations, retaining side-qualified associations, masses and temporal shifts. The common-energy shift is composed with the existing LU affine pullback; its active determinant includes `tau^-3`. At physical adoption, the shared directed ray arithmetic encloses those original equations on the actual canonical and native completed points. Normal displacement and the native host residual each have half of the sampling-accuracy budget relative to the canonical normal radius. Amplitude joints use the same check with zero host residual. Errors trigger physical precision rescue from the retained source, without redrawing maps or partitions. These checks belong to sampling time, including rotations and failed attempts. They do not certify arbitrary CT multipliers, higher raised jets or projected CT-star targets. The #link("../research/advanced_sampling/HOSTED_JOINT_NORMAL_ALIGNMENT.typ")[hosted-joint audit] records the concrete geometry and acceptance gates. The old LMB-specific partition implementation is removed. Explicit `sampling.lmb_basis_ids` overrides preserve exactly the requested basis order; automatic soft-coverage augmentation applies only without that override.

A named standalone cut channel accepts `radial_profile="lu_h"`, or its detailed table form. Warmup fits a normalized log-logistic proposal to the actual runtime LU h-function and compiles its CDF and derivatives alongside that channel\'s proxy program. A positive broad component induces the ordinary raw-radius law independently of the cut root. The map inverts the full mixture CDF, and inverse density uses the supplied raw radius. Named proposals share registered cut geometry without overwriting each other\'s settings. The largest raised order is collected before equivalent cut geometries are deduplicated; this initial fit does not optimize derivative envelopes. Physical h and all residue/CT derivatives remain unchanged. A composed channel may attach this profile to its unique phase-space block. The chart requires the actual LU fixed point to be certified interior; a geometric recentering cannot silently change the meaning of its radial LU scale. The profile is rejected on amplitude channels.

Explicit amplitude E-surfaces populate the production compile context from real catalogue equations, masses and external momenta. Complete native parent LMBs support both full-space and proper-fiber charts. Requested parent order is resolved by reordering a clone of an existing generated basis, preserving its external routing and the catalogue\'s immutable IDs. Active axes follow parent order; an explicit surface shorthand samples its complement first, then the conditional active block. Dependent products are rejected. Ambiguous shifts, unsupported frames and rank-deficient requests fail clearly.

The same implicit map prepares a center and existence status once per forward or foreign inverse, outside root iteration. Two-energy rank-one fibers use an analytic minimum including fixed energies, boosts and unequal/zero masses. General fibers use the existing SOCP center solver only as a candidate source; the actual native energy residual must certify an interior point. Conservative sign tolerances include cancelling routed inputs. Empty-complement preparation is frozen at warmup. Certified absence or explicit pinched status retains a normalized full-support fallback; a failed solve or numerically ambiguous sign requires a typed numerical error. Production preparation now uses one fixed-Arb proposal; only subsequent physical materialization/evaluation participates in precision rescue. Frozen-binding successes and typed numerical failures use the same native runtime caches. Warmup requires the canonical Arb binding for all graphs; optional native component failures cannot reject a valid canonical draw. This is not an interval certificate. Automatic discovery remains unfinished. The combined X2 numerical gates for conditional cuts and frozen-binding rescue pass, including the broader 141-test core/API run. The selected UV-finite two-loop amplitude benchmarks are recorded in #link("../research/advanced_sampling/AMPLITUDE_BENCHMARK_CANDIDATES.typ")[the amplitude study];. The kernel certifies residuals and rejects unrepresentable inverse points without clipping. The shared radial profile focuses on both signed distances to the threshold with a common analytic inverse. Native kernel tests cover Quad and arbitrary precision. Components, affine/composed maps, conditional cut data, eager/dual evaluations and foreign-density partitions now retain native values throughout. The existing eager Jacobian evaluator can select ordered active parameter columns while other prepared parameters remain fixed. Its single compiled first-derivative program receives the statically zero components of the identity seeds through Symbolica\'s `Dualizer`. This prevents an unused derivative of a prepared-only expression such as `sqrt(m)` at `m=0` from contaminating a valid active derivative. Requested singular derivatives and singular intermediates involving active inputs still produce typed numerical errors. The constructor change preserves the stored evaluator layout and the empty-mask behavior of existing physical evaluators. The graph-independent shared-energy joint component supplies a compact three-dimensional full-circle map for two energy sums with one common routed energy. It uses that eager Jacobian owner, a directed MPFR enclosure certificate for the residual disk, and an inverse density evaluated at the supplied point. The existing block resolver now treats two distinct plain surfaces under `intersect` as one three-dimensional block, with qualifiers on the pair. Their separate equations remain in the same geometry key, and the canonical program cache compiles one neutral joint expression for all participating channels. The graph matcher has prior coverage against actual routed amplitude equations. The current source adds an amplitude binder for that supported class, retaining the native shared-energy translation. Empty-prior blocks bind once at warmup through static `Affine`; only nonempty priors use the existing ordered embedding. Its trial disk radius is `e_cm * sampling.b`, with isotropic normal scale one; the component\'s conditional fallback retains the existing normalized radial map. The binder and source transport pass combined numerical gates: all-18-orientation kite Gaussian/moment acceptance, actual physical map/partition consistency and controlled Double-to-Quad body retry, plus the generated triangle\'s empty-prior binder, nonzero offset, compact roundtrip and production reference evaluation. These joint fixtures clear rotation probes; nonidentity-probe behavior and GL638 runtime improvement are not established by this milestone.

The original-source preparation phase now retains the complete 1000-bit Arb draw, including graph-parent coordinates, combined factors and selected host records. Selected, summed and default/LMB routes prepare before norm-based lane selection or physical targets. The existing Gamma/Discrete owners hold completed rows with their event-group boundaries; the old cube-only summed placeholders and physical LMB reinterpretation are removed. Physical attempts materialize directly from the same anchor without maps, partitions or sampling root solves. Direct selected momentum input prepares only its partition; direct unselected input needs no map phase. Norm and debug access do not redraw. `EvaluationMetaData` retains discrete component decisions and root diagnostics, while `EvaluationSource` borrows the immutable draw. Warmed maps have no mutable draw cache. This fixes the continuous-map switching gap described in the #link("../research/advanced_sampling/CANONICAL_DRAW_SOURCE.typ")[canonical source contract];. The correction passes 214 core tests and three API gates, including the complete generated-amplitude and saved-state acceptance runs. Hosted joint binding and completed-point normal accuracy were the next cross-section gate; their later implemented ownership is described above.

Physical cross-section overlaps now also have an immutable source authority. Before lane selection, the existing Gamma row traversal prepares complete accepted-cut representatives from the canonical Arb point through the physical LU root owner, then chooses identity-frame SOCP/forced centers once. Bare raw input retains its real absent channel; selected, summed and nonmaster rows use the same physical convention. Selector-only events carry the source\'s actual optional channel ID, so an explicit channel-ID predicate may change the accepted cut set. Geometry is channel-independent for the same accepted physical cuts. Reference targets and CT-off skip this preparation; amplitudes retain their existing center behavior.

Rows retain the original-frame overlap authority across casts and rotations. Native consumers preserve per-cut complements, LU/alpha solves and raised packets, validating complete group membership, solve signatures and ordered parents before using the retained centers. A mismatch fails without choosing a replacement center. Both consumers promote binary64 center bits exactly before rotating in native arithmetic. Canonical event/runtime and radial history are isolated from physical probes and restored on errors. The `canonical_physical_preparation_time` metadata field exposes this work as a subset of physical time, counted once; map/adoption costs retain their separate sampling attribution. The #link("../research/advanced_sampling/CANONICAL_CT_CENTERS.typ")[center audit] records the generated gates and successful actual GL638 failure replay, without claiming amplitude covariance or numerical coverage of every selector/group.

The geometry and timing source passes 175 selected core tests, including the complete 8192-draw kite rerun and its physical/retry checks, plus three API gates. `Esurface::solve_lu_cut` returns the native ordered energy ray together with its root, and physical raised eta jets evaluate that ray. Explicit LU now routes before scaling; other scalar CT/fiber evaluators retain their prior operation order and seed convention. Those gates predate the native host handoff described below and do not validate its adoption checks. Three #link("../research/advanced_sampling/GL638_LU_RAY_REPLAY.typ")[GL638 LU-ray replay evaluations] retain all 936 orientations, six physical event weights and the unchanged saved state. They establish a bounded physical regression, not improved precision or variance; GL638\'s six simple cuts do not exercise higher raised eta derivatives.

The accompanying timing change passes core/API checks. `parameterization_time` accumulates map/source/partition and replay costs across attempts; canonical preparation is charged once. `integrand_evaluation_time` counts actual target bodies across all lanes and probes, excluding maps and nested host-adoption verification charged to sampling time. Evaluator and event timings are subsets removed from the benchmark\'s residual integrand row. Final metadata retains accumulated costs, including replays after its earlier result snapshot. These are measurement counters, not evidence that the GL638 10% budget is met. The inclusive evaluator timing update counts every physical evaluator wrapper call across probes and rescues, including completed calls before a failure and the raised-threshold IFT alpha helper. Earlier reports counted only the primary rotation at the first precision and cannot supply this new E split. The `canonical_sampling_preparation_time` subset uses the same elapsed interval as its contribution to S; the existing canonical physical preparation remains a subset of P. Bench rows subtract these subsets once and reject inconsistent subset durations. These additions are per-evaluation metadata only; the serialized `StatisticsCounter` and integration checkpoint layouts are unchanged. Multiworker benchmarks must report summed worker costs separately from elapsed wall time.

The native host handoff is implemented and its combined gates pass. The existing runtime context owns one native record vector; the initial selected forward or direct inverse freezes its authoritative prefix before partition inverses. Only identical plans and exact represented priors reuse a record. Selected/summed samples lend that prefix to the physical owner. The subsequent source correction materializes it directly from Arb, and rotations retain the original-frame payload. Metadata contains discrete policies and root history, with canonical preparation history isolated. Graph-group masters authenticate their records before ordinary foreign graphs receive an empty payload. Physical adoption rotates the ray once and authenticates its parent/host against the selected immutable channel, then compares it with independently routed completed data using a fixed directed 2048-bit scalar/root/first-derivative/inverse-volume check. This does not certify vector-coordinate error, higher raised jets or full conditional density. The final source passes 184 unique core tests and three API gates (106.024 s). The conditional fixture\'s nonidentity physical totals agree within the fixed 1e-8 tolerance; probes retain no event buffer/counts under existing ownership. The complete 8192-draw kite passes in 1504.601 s with its 6%/8% reference bounds; this unoptimized complete-fixture timing is not a production budget measurement. Combined checking (12.56 s), clippy (57.76 s; no changed-line warnings), formatting and diff checks pass. The optimized GL638 10% budget remains unmeasured.

That native-host milestone preceded hosted joint binding. The subsequent fixed-external routing, retained canonical law and original-normal adoption checks are described above; CT-star maps still need to consume the actual common center and projection alpha through their own dependency-certified handoff. Map coverage is `Full` or `Restricted`, independently of `requires_context` and the new `requires_proposal_policy` preparation flag. Certified inverse exclusion propagates through the existing compositions, embeddings and affine frames as `None`, becoming zero in the existing partition. Numerical uncertainty and underflow remain errors. Exact and proxy partitions gate restricted channels by their actual support, and bridge admission requires an explicitly selected full-support sibling. A direct-momentum request outside its selected channel gives an explicit diagnostic; generated selected points must still have positive density and satisfy forward/inverse consistency. Typed numerical errors distinguish unrepresentable derived coordinates/scores from invalid original inputs. The implicit root uses the existing native safeguarded solver, with callback errors preserved. Production retries materialize the retained canonical draw at the requested precision with matching improved native external data, reusing its mapped points and foreign-density partition. A typed physical failure discards partial values and events before retry. Precise output retains native values; ordinary output rejects unrepresentable finite values instead of silently zeroing them. Root uncertainty relative to the sampled threshold distance and density accuracy still needs its stronger certificate. Derived-expression mass evaluation, Gaussian-body underflow and final outer-grid range handling remain separate precision limits. The bridge checks the selected forward determinant against its inverse density at the actual mapped point. It reuses a selected score only when that channel actually uses its map density; OSE and user-proxy overrides evaluate the selected inverse independently. Warmup assigns one tenth of the strictest matching stability-level relative tolerance as its density budget; standalone constructors use native square-root epsilon. Failure is a typed numerical error. This is numerical proposal consistency, not a rigorous enclosure of the physical surface or of floating-point errors.

Channel weights are resolved per canonical entry: the global `sampling_channel_weight` is overridden by an optional named `channel_weight`. Genuine OSE weights are available for complete ordinary LMB maps, using raw master-frame edge energies and the dimensionally consistent score `E_cm^(-3L) product(E_cm/E_e)^alpha`; surface maps retain their exact densities unless an explicit compatible user proxy is selected. The existing positive-score evaluator and support-aware partition own all strategies. `SamplingChannelSelection` owns `weight` and `alpha` together; no separate legacy LMB-weight owner or enumeration is retained. See the #link("../research/advanced_sampling/OSE_CHANNEL_WEIGHTS.typ")[OSE settings and acceptance contract];.

Physical and Gaussian-reference targets share graph traversal, canonical sampling and native physical precision retries from retained draws. Reference values and raw-frame moments each pass the existing scalar stability check before joint reporting; unresolved or unrepresentable results fail acceptance. Summed reference evaluation aggregates values and raw-frame moments across mapped points before statistical squaring. Graph/orientation probabilities and outer map/grid factors are applied once; invalid draws fail acceptance. Explicit discrete coordinate selections expose partial partition contributions, which require summation or inverse-selection probabilities for unit normalization. The ordinary `integrate` workflow also accepts a `reference_gaussian` descriptor. CLI syntax is `--reference-gaussian '{"width":300}'`; Python callers pass the same command through `api.run(...)`. The retained `reference_gaussian=(300.0, [])` tuple helper, whose second item is an optional shift, belongs to the disabled typed integration entry point. Its existing complex accumulators report Gaussian normalization in Re and `|K|²` times the Gaussian, divided by the known raw second moment, in Im. Both expectations are one. This is acceptance output, not a physical imaginary part; the raw moment scale is reported separately. An omitted center is resolved to zero in the loaded process\'s spatial dimension. Physical selectors and observables must be disabled for this substitution. The initial binding requires one common nonzero loop count across each slot\'s loaded graphs; it does not infer a smaller dimensional scope from graph filters.

The reference stays on the existing integration slot and is passed explicitly to worker evaluations. The same native stability loop checks its value and raw moment before reporting. No second integrand, grid or statistics engine is introduced. Workspace manifest version 3 records each resolved descriptor and the observable convention; resume requires the same reference request and checks it before restoring runtime settings. Older integration workspaces require `--restart`. Generated `gammaloop_state/` data and the binary integration checkpoint layout are unchanged by this reporting extension. See #link("../../ADVANCED_SAMPLING_PLAN.typ")[the implementation plan] for the remaining harness and channel milestones, and #link("../research/advanced_sampling/LU_H_MATCHED_SAMPLING.typ")[the LU localization study] for the implemented h-matched profile and its remaining optimization limits.

=== 3.2 Differential event model
<32-differential-event-model>
Events are stored in precision-homogeneous containers:

- `GenericEvent<T>`
- `GenericEventGroup<T>`
- `GenericEventGroupList<T>`

This avoids enum-based mixed-precision layouts that would otherwise
force storage to be sized by the largest precision variant.

Each event carries:

- kinematics
- cut metadata
- one fully normalized event `weight`
- `additional_weights`

Observable-specific entry reweighting remains internal to the observable
runtime; it is not stored on the event.

`additional_weights.weights` is a generic `BTreeMap` keyed by lightweight
identifiers such as:

- `FullMultiplicativeFactor`
- `Original`
- `ThresholdCounterterm { subset_index }`
- `AmplitudeThresholdCounterterm { esurface_id, overlap_group }`
- `AmplitudeThresholdCountertermVariant { variant_id, esurface_id, overlap_group }`

Counterterm weights are stored with the sign with which they contribute
to the final event weight. Cross-section events use
`ThresholdCounterterm`; amplitude events use
`AmplitudeThresholdCounterterm`. The fully normalized event weight is
therefore reconstructed as:

`(Original + sum(the event's counterterm entries)) * FullMultiplicativeFactor`.

This is populated only when
`settings.general.store_additional_weights_in_event = true`.

Directive-enabled events can also retain an addable `threshold_counterterms`
decomposition. A component ID refers to the graph's static registry, while an
occurrence identifies its amplitude surface/overlap or local-unitarity overlaps
and raised orders. Each component records the multiplier values, effective
multiplier, fully normalized weighted contribution, and bare contribution with
only the user multiplier omitted. When an exact-zero multiplier skips the bare
evaluation, `bare` is absent, `weighted` is zero, and `evaluation_skipped` is true.
The Rust and Python event surfaces expose the same component and provenance
boundaries. JSON omits absent decompositions and bare weights and retains the
occurrence's `kind` tag. Positional binary serialization keeps every optional-field
marker and uses enum discriminants, preserving round trips through both formats.

For physical W/Z final states, generation includes the model-declared covariant vector, Goldstone and ghost cut states. Event construction preserves their momenta and maps their PDGs to the requested physical vector before selectors or observables run. The same event receives its bare, UV and threshold weights, so measurement functions respect the complete gauge sum. The representative map is persisted with each cross-section graph term. Ambiguous requests mixing an explicit unphysical state with its physical-vector sector are rejected; separate diagnostic requests retain their original PDGs. See the #link("sm-conventions-audit.typ#electroweak-virtual-and-cut-state-gauge-contract")[electroweak gauge contract] for the required graph completeness and treatment of intentional subsets.

Serialized event `additional_weights.weights` is a list of `[key, value]` pairs. This preserves structured counterterm identifiers in JSON, including `AmplitudeThresholdCounterterm { esurface_id, overlap_group }`, and supports round trips without converting identifiers to display strings. `inspect --json-output` retains the same event weights and evaluation as ordinary inspection. Consumers of event JSON must read the pair list, including `[]` for no weights.

=== 3.3 Event generation policy
<33-event-generation-policy>
Event generation is controlled by three related runtime decisions:

- `should_generate_events()`
  - true when `general.generate_events = true`
  - or active selectors exist
  - or observables exist
- `should_buffer_generated_events()`
  - true when `general.generate_events = true`
  - or observables exist
- `should_return_generated_events()`
  - true only when `general.generate_events = true`

This means:

- `general.generate_events = false`
  - selector-only evaluations still build temporary events for selector
    checks
  - selector+observable evaluations still build temporary events and
    buffer accepted identity-rotation events long enough to fill
    observables
  - those temporary events are cleared from returned sample results and
    integration event outputs afterwards
- `general.generate_events = true`
  - events are built, buffered, and returned even when no selector or
    observable is configured

The standard Rust/Python sample-evaluation helpers may temporarily
override the returned-event behavior per call, but internal event
generation still follows this same policy after that override is
applied.

Accepted LU events are retained as grouped event lists:

- one `EventGroup` per graph-group evaluation
- all accepted cuts from graphs sharing the same `group_id` are merged
  into the same retained event group
- each event stores its concrete `graph_id` and `cut_id`
- incoming PDGs are sourced from the initial-state cut edges

This preserves the correlation structure needed by downstream
observables.

=== 3.4 Observable architecture
<34-observable-architecture>
Observable definitions are shared between selector and histogram
use-cases:

- one quantity definition extracts zero or more `ObservableEntry<T>`
  values from an event
- selectors consume those entries on a single event
- histogram observables consume grouped events and fill per-group
  contributions

Histogram accumulation is based on:

- a histogram-level `sample_count`
- sparse per-bin sufficient statistics:
  - `entry_count`
  - `sum_weights`
  - `sum_weights_squared`
  - `mitigated_fill_count`
- `HistogramAccumulatorState` / `ObservableAccumulatorBundle` for
  mergeable runtime state
- `ObservableSnapshotBundle` for JSON / HwU output and API responses

Each Monte Carlo sample contributes exactly one statistical sample to a
histogram. Contributions from all retained event groups are summed first
per bin, and untouched bins are handled implicitly through the
histogram-level `sample_count` rather than explicit zero fills.

Histogram snapshots therefore remain fully mergeable and
re-constructible into live accumulator state. The Rust-facing histogram
API exposes `merge(...)`, `merge_in_place(...)`, and `rebin(...)`.

Current built-in quantities include:

- `particle` with `computation = scalar | count | pair`
- `jet` with `computation = scalar | count | pair`
- `afb`
- `integral`
- `graph_id`
- `graph_group_id`
- `orientation_id`
- `sampling_channel_id`

Current pair quantities include `DeltaR`. Current scalar projections
include `E`, `CosTheta`, `PT`, `y`, `eta`, `Px`, `Py`, `Pz`, and `Mass`.

Object and pair quantities are ordered before `entry_selection` is
applied. The ordering contract is configurable per quantity through:

- `ordering = PT | Energy | AbsRapidity | Quantity`
- `order = Ascending | Descending`

Family defaults are:

- particle scalar quantities: `ordering = Quantity`
- jet scalar quantities: `ordering = PT`
- pair quantities: `ordering = Quantity`

`leading_only` and `nth_only` therefore operate on this explicit ordered
list, not on incidental source/event enumeration order.

Jet clustering is implemented natively in Rust with IRC-safe algorithms:

- `kt`
- `cambridge_aachen`
- `anti_kt`

and is validated against `fjcore` in tests.

== Persistence and Runtime Artifacts
<persistence-and-runtime-artifacts>
Primary persisted state lives under `gammaloop_state/` (default):

- `state_manifest.toml` (state schema/version marker)
- `model.json`, `model_parameters.json`
- `symbolica_state.bin`
- `processes/` (amplitudes/cross\_sections + integrands)
- `run.toml`
- `default_runtime_settings.toml`
- `global_settings.toml`
- `logs/`

For local experimentation, prefer an isolated path such as
`.local/scratch/<run>/gammaloop_state` to keep repository root output
minimal.

The persistence model is file-system based and intentionally
human-editable for settings/run cards, mixed with binary artifacts for
performance-heavy data.

=== Persistence Compatibility Contract
<persistence-compatibility-contract>
- State format is versioned with `state_manifest.toml` (`version = 9` currently).
- Version 9 combines Symbolica 3 evaluator/rational payloads with native sampling and threshold metadata. It rejects both upstream version 7 and sampling-branch version 8 states; regenerate them rather than relabeling or migrating positional bincode data. Standalone amplitude and cross-section payload versions are 10 and 13, respectively. Version 7 stored exact CFF coefficients as native rationals. Version 6 removed obsolete deferred-integrand fields; version 5 added component-local generated-CFF ownership and prefactor metadata; version 4 added the typed global-prefactor sign.
- State loading and direct overwrite both require exactly the current manifest version; older states must be regenerated, and states from newer binaries require a newer GammaLoop binary.
- A missing manifest denotes an unmanifested folder rather than a legacy state and is never loaded as saved state.
- Process settings history now uses `settings_history.toml`
  consistently; loader still accepts legacy `settings_history.yaml` for
  backward compatibility and migration.

=== Integration workspaces
<integration-workspaces>
Integration resume state is persisted in a dedicated integration
workspace. In a normal writable session its default is
`<state.folder>/integration_workspace` (therefore normally inside
`gammaloop_state/`); a read-only session defaults to
`./integration_workspace[_<state-name>]`, and `--workspace-path`
overrides the default. `manifest.json` records the selected
process/integrand slots, targets, effective model parameters, integrand
fingerprints, training and settings slots, and sampling-correlation
mode. Authoritative completed-iteration state lives in
`state/integration_state.bin`; per-slot settings live under
`integrands/<process>@<integrand>/settings.toml`, and observable resume
snapshots live under `state/observables/<process>@<integrand>/`.
User-facing `integration_result.json` and observable files are derived
snapshots, not the resume authority.
Workspace manifest version 3 also retains the resolved Gaussian-reference request
and the convention for physical absolute contributions at distinct channel points.
Earlier checkpoints require a fresh run; resume verifies the reference request
before restoring settings or completed-iteration state.

== Configuration Architecture
<configuration-architecture>
Configuration is split into:

- Global settings (`GlobalSettings`): generation, logging directives,
  parallelism.
- Runtime settings (`RuntimeSettings`): kinematics, integrator behavior,
  sampling, stability, subtraction, and differential event-processing
  configuration.

Runtime updates normalize the input aliases `lmb_multichanneling`, `lmb_channels`, and `lmb_channel_weight` to `sampling_multichanneling`, `sampling_channels`, and `sampling_channel_weight` before merging with stored settings. These aliases work with files, TOML strings, and key/value updates; one update cannot specify both spellings of the same field.

Differential runtime configuration currently includes:

- `general.generate_events`
- `general.store_additional_weights_in_event`
- named `quantities`
- named `selectors`
- named `observables`
- `integrator.observables_output`

This split is good: generation-time and evaluation-time concerns are
separated and independently serializable.

The maintained #link("../../../products/gammaloop/latest/guides/diagnostics/")[logging and
diagnostics guide] documents the user-facing startup precedence, filter language, tag
vocabulary, and sink-routing contract implemented by these settings.

== Concurrency Model
<concurrency-model>
Concurrency is explicit and use-case scoped:

- Generation and compile thread pools use configurable thread counts.
- Integrator parallelism is controlled via runtime/global settings.
- Some loops over processes/integrands remain sequential at
  orchestration level while heavy operations inside are parallelized.

For the differential pipeline:

- each worker/integrand clone owns its own mutable observable
  accumulator state
- local worker results are merged back into the master integrand after
  chunk evaluation
- distributed batch mode can return either grouped events or
  serializable observable accumulator bundles
- observable snapshots are written only from the merged master state

== Testing Architecture
<testing-architecture>
Testing is organized into layers:

- Unit tests in modules across `crates/gammalooprs/src/`.
- Integration test crate: `tests/` with reusable harness
  (`tests/src/lib.rs`).
- Additional Rust integration tests in `crates/gammalooprs/tests/`.
- Snapshot-heavy coverage for symbolic outputs and numerics.

Differential coverage is currently split into:

- `tests/tests/test_clustering.rs`
  - compares the native jet clustering implementation against `fjcore`
- `tests/tests/test_differential.rs`
  - validates grouped event propagation
  - selector behavior
  - observable filling and merge semantics
  - JSON / HwU snapshot output
  - API-facing result shaping using handcrafted fixtures
- `tests/tests/test_evaluation_api.rs`
  - real generated-process end-to-end coverage for Rust
    `evaluate_sample(s)`
  - graph/cut grouping and incoming-PDG propagation checks
  - selectors / observables / event-group retention
  - batch-global observable snapshots for `evaluate_samples`
  - per-sample evaluation metadata
  - x-space vs momentum-space Jacobian behavior
  - minimal-output mode
- `tests/tests/test_runs.rs`
  - real generated-process end-to-end coverage for `integrate`
  - integration-workspace observable output in JSON / HwU formats
  - repeated `save dot` overwrite behavior
- `tests/tests/test_python_api.rs`
  - subprocess-based end-to-end coverage for the public Python API
  - returned event grouping / ids / incoming PDGs
  - selectors / observables / returned metadata
  - minimal-output mode
  - x-space batch and momentum-space sample coverage
  - requires a separately prepared Python environment where
    `import gammaloop` already works

== Architectural Strengths
<architectural-strengths>
- Clear separation between application shell (`gammaloop-api`) and
  domain core (`gammalooprs`).
- Rich, serializable settings model with schema generation hooks.
- Strong process abstraction (`Process`, `ProcessCollection`,
  `ProcessIntegrand`) that supports amplitude and cross-section
  workflows.
- Robust stability pipeline with multi-precision escalation and rotation
  checks.
- Differential event-processing is now integrated into the core
  evaluation path without compromising the stable-weight selection
  logic.
- Observable accumulation is mergeable across rayon workers and
  distributed batch outputs.

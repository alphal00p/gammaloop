# gammaLoop Current Architecture

## Scope
This document describes the current, implemented architecture of this repository.

The workspace is organized around two main Rust crates:
- `gammalooprs` (`crates/gammalooprs`): core physics/domain logic, graph processing, integrand construction, evaluation, and integration.
- `gammaloop-api` (`crates/gammaloop-api`): CLI, REPL, Python bindings, command parsing, and persisted state orchestration.

## High-Level Architecture
At a high level, gammaLoop uses a layered architecture with a stateful application shell.

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

## Main Components

### 1. Entry Points and Interfaces
- CLI entry point: `crates/gammaloop-api/src/cli.rs`.
- Main app orchestration: `crates/gammaloop-api/src/lib.rs` (`OneShot::load`, `OneShot::run`).
- Python module entry: `crates/gammaloop-api/src/python.rs` (`#[pymodule(name = "_gammaloop")]`).
- Python package shim: `crates/gammaloop-api/python/gammaloop/__init__.py`.

### 1.1 CLI architecture
- The CLI is state-centric: every session operates on one active state folder and one active in-memory `RunHistory`.
- Startup is orchestrated by `OneShot`:
  - resolve the state folder from CLI / boot card / defaults
  - load or create the state
  - load persisted settings and run history
  - optionally apply a boot card before entering the REPL or running a one-shot subcommand
- Command execution is centralized in `crates/gammaloop-api/src/session.rs` (`CliSession`):
  - top-level commands
  - boot-card replay
  - nested `run` execution
  - command-history recording policy
  - command-block-definition mode
- The REPL layer in `crates/gammaloop-api/src/repl.rs` is state-aware:
  - prompt text reflects the active state settings
  - tab completion is driven from the active state, active command blocks, and clap argument metadata

### 2. Application State and Command Model
- Central mutable app state: `crates/gammaloop-api/src/state.rs` (`State`).
- Command history and run cards: `CommandHistory`, `RunHistory`.
- Command dispatch: `crates/gammaloop-api/src/commands/mod.rs` with concrete subcommands in `commands/*`.

The command model is stateful by design: commands mutate a long-lived `State` that can be saved and resumed.

### 3. Domain Core (gammalooprs)
- Root module wiring: `crates/gammalooprs/src/lib.rs`.
- Model and parameters: `crates/gammalooprs/src/model/mod.rs`.
- Graph domain: `crates/gammalooprs/src/graph/mod.rs` and submodules.
- Process orchestration: `crates/gammalooprs/src/processes/mod.rs`, `crates/gammalooprs/src/processes/process.rs`.
- Amplitude and cross-section pipelines:
  - `crates/gammalooprs/src/processes/amplitude.rs`
  - `crates/gammalooprs/src/processes/cross_section.rs`
- Integrand abstraction and implementations:
  - `crates/gammalooprs/src/integrands/mod.rs`
  - `crates/gammalooprs/src/integrands/process/mod.rs`
  - `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`
  - `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs`
- Integration engine: `crates/gammalooprs/src/integrate/mod.rs`.
- Global/runtime settings: `crates/gammalooprs/src/settings/mod.rs`, `crates/gammalooprs/src/settings/global.rs`, `crates/gammalooprs/src/settings/runtime.rs`.
- Differential event/observable pipeline:
  - `crates/gammalooprs/src/observables/events.rs`
  - `crates/gammalooprs/src/observables/observables.rs`
  - `crates/gammalooprs/src/observables/clustering/*`
  - `crates/gammalooprs/src/integrands/evaluation.rs`

### 3.1 UV renormalization

UV generation defaults to the disconnected-capable hedge-poset backend and
also supports a legacy DAG-forest backend and a comparison mode. Scheme policy
remains on each Spinney while compute nodes store typed local, integrated, and
cut-dependent results. Four-dimensional disconnected counterterms factorize
over complete component counterterms; three-dimensional disconnected terms
replay component-local operations from their common root so CFF structure is not
multiplied as though it were scalar. The maintained sign, projection, marker,
and backend-boundary invariants are documented in
[`uv-renormalization.md`](uv-renormalization.md).

### 3.2 CFF production and numerator-energy ownership

GammaLoop owns production graph/source construction, UV orchestration, exact
source mapping, and evaluator preparation. The `three-dimensional-reps` crate
owns the shared CFF algebra. The `3Drep` command and feature-gated eager
evaluator are diagnostic tools, not production contracts: GammaLoop may prepare
their inputs, factors, and expressions differently.

All CFF power and capacity questions, including numerator and finite-pole
powers of repeated occurrences, are expressed solely in physical EMR/source-edge
energies. LMB coordinates describe momentum routing and are never consulted to
identify, cap, or substitute an energy power. Production capacity analysis
accepts physical `Q(edge, index)` atoms and rejects `K(loop, index)` until its
producer normalizes it with physical edge provenance.

Both direct local-3D modes first perform the complete loop-energy integration
and build the complete/global CFF expression; the UV Taylor operators then act
on that CFF expression. They use the same Taylor-transformed CFF bodies.
Writing the generalized residue map as `{ k -> C_k }`,
`explicit_orientation_sum_only=false` is `sum_k sigma(k) C_k`, with `sigma(k)`
the one-hot selector for the complete residue-map key. The Taylor operator is
applied independently to every keyed body and leaves that selector opaque.
`explicit_orientation_sum_only=true` only replaces each selector by one and
explicitly sums the same bodies. Neither direct mode reconstructs
or projects completed local-4D Taylor structures.

Only the projected local-4D route represents completed terms with raised
propagators by source-backed occurrence graphs. The original `EdgeIndex`
survives the Taylor operator in the typed denominator wrapper. Disjoint-set
contraction of absent edges in the original graph then constructs the cograph
and UV source minors, and every
occurrence inherits the endpoints of its `source_edge` in the appropriate
minor. Only repeated occurrences of that same source edge subdivide its
incidence into a serial dotted chain. Exact momentum signatures are normalized
up to sign to validate denominator equality and record repeated occurrences;
they never infer endpoints or merge physical owners. The rank solve at this
source-topology boundary selects a unique `+/-` routing sign modulo the opposite
source domain on those
already fixed endpoints. There is no incidence/Kirchhoff reconstruction and no
external-balance synthesis. The raw `+/-Q` sign remains available to the
numerator mapper, while a post-construction Graphica pass canonically relabels
nodes and exact edges for deterministic cache keys.

Projection prepares each genuinely outer additive Taylor term once, then
compares at most three certified assignments by actual native source-map count.
Generated expressions are reused only for equal canonical topology and
identical per-occurrence capacity. Every term keeps its selected factorized
assignment and matching generated payload.
Equal physical energies do not authorize the cache to redistribute bounds
between their distinct occurrences. Independent requests never contribute a
shared capacity maximum, and no preliminary registration pass is needed.
The real degree-one triangle regression retains the natural UV-owner
multiplicities `(1,1,1)`, `(2,1,1)` and `(1,1,2)` as separate denominator
topologies. Terms with the same denominators may share a factorized additive
numerator; collecting a common denominator across different topologies would
manufacture positive denominator factors and inflate the required CFF rank.
Genuine positive typed denominator factors already present in a numerator
remain supported, with their original ownership and CFF lower sectors intact.
Non-vacuum exact sources retain pure-external boundaries as explicit
source-crown hedges. Future on-shell two-point insertions such as `(m,0,0,0)`
require an explicit fixed-boundary payload, not topology reconstruction.
This exact-source reconstruction and bounded dispatch search is exclusive to
the projected local-4D route; direct local-3D and LTD dispatch are unchanged.
Original factors retain their own occurrences; newly denominator-derived hard
factors may use only their line's serial copies. The descending rank envelope
orders proposals, preferring `(4,2,2)` to `(4,3,1)` while the original quartic
stays fixed. The hard search compares that baseline with at most two
deterministic alternatives, each changing a single owner's assignment. New soft
Taylor factors retain explicit crown provenance until outer assembly, where the
completed child coefficient is multiplied by the untouched outer numerator.
The soft search proposes at
most three certified off-shell cograph routings, including every fixed external
shift.

Both searches select the smallest actual native generated source-map row count
(`generated.expression.orientations.len()`), with proposal order breaking ties.
This counts stored operational rows after the core's own emission/coalescing,
including retained zero rows, before surface conversion and cut/host selection.
It is not a degree-product estimate or a global union of keys across Taylor
terms. Final selected host/map branches can have a different count. Existing
rank/Pareto and basis restrictions define a bounded proposal class, so this is
best-of-three selection, not a global optimum or a bound on generation time/RAM.

The hard cache admits full payloads only for winning keys and retains small
count-only records for contenders under the same complete source/options/bounds
key. The soft path has no persistent count memo: it holds an incumbent and
challenger, drops the loser, then performs surface conversion, cut/hosting,
reporting and numerator mapping only for the winner. Existing root-expression
reuse bypasses this search instead of regenerating its established capacity.
The chosen factorized expression and immutable plan supply both the generated
capacity and numerator mapping. No numerator expansion, on-shell routing
identity, or change of original ownership is involved.
The full design, sign argument for `D(Q)=D(-Q)`, and concrete production fixtures
are documented in
[`exact-powered-denominator-cff-lifting.md`](exact-powered-denominator-cff-lifting.md).

The production numerator remains factorized. Degree analysis traverses its
factors without expanding them, each UV step attaches only newly owned factors,
and final assembly attaches outside and global factors exactly once. When one
outer CFF serves independently evaluated residue branches, its capacity is the
per-edge maximum of their separately analyzed factorized products, including
the common remaining numerator. Equal branch bodies need only one rank analysis;
opposite bodies remain independent and cannot cancel each other's capacity.
The analysis streams degree maps without constructing a tagged symbolic sum.
Initial-cut and tree-edge energies remain outside this internal CFF capacity
domain; stored production-root CFF reuse retains its existing boundary. For higher
power projection, the term parser splits only the completed expression's outer
Taylor sum when its addends carry separate denominator topologies. Nested
numerator sums and positive typed `GS.den` factors remain factorized and are
never cancelled against denominator occurrences upstream; generalized CFF owns
the resulting pinches and lower sectors. For higher powers,
interpolation may replace an EMR energy by `a*M`, where `a` is a signed
integer and `M` is the common auxiliary CFF numerator-sampling scale. This is an
EMR substitution, never an LMB rewrite. Production evaluator parameter lists
always include `M`, and runtime settings require `M != 0` for every process.
The physical result is invariant under changing its nonzero value.

The signed sampling coefficient is neither a physical edge-direction sign nor
the runtime residue-map-key selector `sigma(map_id)`. Pole signs enter the generation-time affine map, after which an
orientation entry stores `a` in `LinearEnergyExpr::uniform_scale_coeff` and its
interpolation weight stores the compensating inverse power in
`CFFVariant::uniform_scale_power`. Multiple numerator maps may therefore share
the same coarse orientation. Their complete loop/edge energy maps and distinct
`numerator_map_index` values remain attached to the factorized numerator. Each
such entry has its own complete map key and therefore its own `sigma(map_id)`;
entries must never be merged merely because their physical edge directions
agree.

Integrated finite UV terms retain their exact source-local EMR maps. Production
map-key selectors partition complete generalized residue-map entries but do not
replace those maps.
In orientation-local direct-3D generation, shrinking a UV subgraph preserves
the orientation selected on every surviving outer edge. Among the cut-valid
full orientations which extend that outer assignment, one deterministic inner
orientation is selected to host the complete integrated finite counterterm; the
other inner extensions receive none, so summing orientations counts the addback
exactly once. In explicit-sum direct-3D generation the same reduced residue is
kept once without a selector. Projected local-4D counterterms are likewise
selector-free: each completed exact source owns its full source-local
orientation sum, independently of production-orientation IDs.
After that sum, a projected sector stores only its coefficient and frozen
localizing factor. Its child map lengths need not agree across Taylor
topologies. Actual outer branches keep their own maps until the coefficient
and remaining numerator have been mapped, then sum by cut order. Typed zeros
retain the allowed cut keys without requiring a map; direct and integrated
localization retain their own production-host rules.
The empty UV forest is the ordinary factorized production root in both local-UV
routes; the expanded-4D setting changes only proper, nonempty UV nodes.

Final assembly uses one `Integrands` map of cut-indexed factorized expressions.
UV markers and tensor replacements act on those expressions before the ordinary
evaluator preprocessing; there is no parallel deferred-body state. Direct-3D
sequential forest construction composes the same local and integrated replay
operations used by disconnected forests. Integrated localization registers its
source surfaces before either Taylor branch reads the graph, and its normalized
localizing factor remains inert under subsequent Taylor operations.

The shared CFF core also returns its connected-loop and pure
duplicate-denominator global sign as typed metadata. GammaLoop consumes that
bridge exactly once for root, reduced, and exact production CFF sources,
cancelling the shared-core-local uniform convention before the physical
Minkowski measure `i^L/(2*pi)^(3L)` is applied. Integrated UV addbacks use the
same convention, with additional Vakint normalization `1`.
The NLO acceptance layer independently generates orientation-local direct 3D,
explicit-sum direct 3D, and projected local 4D with local and integrated UV and
threshold counterterms. It compares complete GL0/GL2 values at a common native-
Arb point in all three routes; the fast Monte Carlo tests integrate explicit-sum
3D. DD acceptance checks the inclusive `(alpha_s/pi) * LO` correction,
graphwise UV-mass and localization-scale independence, cancellation of total
renormalization-scale dependence, opposite GL0/GL2 squared-scale logarithms,
and the physical and projected EMR energy bounds. TT acceptance uses the
fully-MSbar scheme, without on-shell counterterms.

Direct-photon benchmarks use the off-shell spin projector `-g^(mu nu)` and no
picobarn conversion. The inclusive lepton-process targets instead use the
Eq. (7.1) normalization `2(4 pi alpha)/(3 Ecm^3)` and the conversion to picobarns;
individual lepton-process graph components are not assigned the unconverted
published photon targets. Acceptances require positive real LO and the signed
real NLO correction, checking the imaginary component separately. Diagram
contributions can have either real sign. Current validation results and
measured timings belong in the accompanying test evidence and PR.

The raw forward graph already contains inverse-process UFO vertices with their
Hermitian-partner spin structures and couplings. An additional RHS numerator
adjoint would conjugate those structures a second time. Marking vertices and
virtual propagators, together with the reversed RHS virtual contours, instead
gives `(-1)^C_R`, with boundary hairs retained in the RHS component count. The
complete LU residue factor is `2*pi*i*(-1)^C_R`, hence `-2*pi*i` for connected
RHS, while retaining the cut propagator numerators exactly once. Bare, local
and integrated UV, and threshold branches use this same cut-group factor.
Runtime subtracts threshold helpers with left `-i*pi` and right `+i*pi`
integrated coefficients; iterated terms use their product. See
[phase conventions and the independent conjugation audit](phase-conventions.md)
for the cut-line `i` cancellation, intrinsic complex CKM phases, and the
unsupported complex-mass-scheme boundary. There are no conjugation modes.

The scalar local-equivalence matrix is generated from the scalar model rather
than from hand-built graph data. Its lanes without an additional numerator
retain the full UFO Feynman rules; companion probes use only Feynman-rule-local edge factors, and
there is no graph-specific production branch. The matrix enables local UV,
integrated UV and threshold counterterms while comparing all three local-UV
routes, including native-Arb checks. `just test_LU_scalar_xs` includes the slow
cases, uses release compilation by default and stops on the first failure.
Near-zero route comparisons retain the documented f64-input floor alongside
native-Arb comparisons; that input floor does not establish arbitrary-precision
agreement. The curated suite includes the six base scalar graphs GL00, GL02,
GL04, GL08, GL09 and GL24. Each profiles direct local3D separately for every
complete residue-map key and projected local4D after the complete residue sum.

UV profiling defaults to `only-divergent`: every expected cycle union with
DOD >= 0 is tested using the generation LMB when suitable, otherwise the first
suitable basis in the deterministically sorted complete LMB list. The
exhaustive `all` mode is opt-in. Amplitude and LU inputs share this behavior,
with graph and Cutkosky-cut selectors for LU profiling and a colored final
failure summary. Per-key profiling is defined for orientation-parametric,
localized direct local3D. Selector-free explicit-sum direct local3D and
projected local4D are summed representations and reject that request.

### 3.3 Tensor-network contraction order

Evaluator construction parses the factorized numerator into a Spenso network,
aliases large scalar references, and contracts its tensor products before
resolving those aliases. The default `intermediate_cost` preset selects
`MinIntermediateCost`, a configuration of the existing `ContractionStrategy`
implementation. It changes pair selection, preserving the scalar, trace,
library, disconnected-product and final cleanup rules.

The greedy score first estimates copied symbolic payload, then symbolic
normalization work, before using sparse support, tensor volume and rank as
later criteria. For a candidate pair, let `J` count the estimated nonzero entry
products with matching contracted coordinates. Let `b` and `t` be each
operand's mean entry bytes and top-level terms, rounded up. The two leading
scores are `J * (b_left + b_right)` and
`J * (b_left * t_right + b_right * t_left)`. Sparse support supplies the join
count without multiplying coefficients. Otherwise the fallback uses the
Cartesian entry count, capped by output and contracted-coordinate volume for
single tensors. That cap is not applied to lazy sums, whose terms may share
coordinates. The existing bounded output-support join and deterministic
operand-order ties remain in use.

These scores estimate work before cancellations. They are not physical memory
bounds or a prediction of the globally best contraction sequence. They use
existing structure and scalar-size profiles, with saturating arithmetic; the
planner performs no speculative contractions or numerator expansion. A higher
rank intermediate can therefore be preferred when its estimated cost is lower.
Performance claims require identical-input measurements of both time and memory.

The previous sparse-aware preset and all other presets remain explicit choices.
To restore the preceding production configuration in a run card, use:

```toml
[cli_settings.global.generation.evaluator]
spenso_execution_mode = ["Sequential", "MinResultRank"]
tensor_network_contraction_order = "sparse_atom_aware"
```

Within a running session, switch the preset with the existing command
`set global kv global.generation.evaluator.tensor_network_contraction_order=sparse_atom_aware`.
Use `intermediate_cost` to select the default again. The preset applies to
`ContractionMode::MinResultRank` with the supported sequential executors;
explicit `SmallestDegree` selection is unchanged. The existing
`bnl_evaluator_atom_mwe` diagnostic exposes the same choices through
`--contraction-order`, using hyphenated names such as `intermediate-cost`.

## Lifecycle and Data Flow

### 1. Startup
1. CLI parses command input (`OneShot::parse_env_with_capture`).
2. Initialization runs (`initialise()`), installs hooks and symbol registries.
3. State is loaded from state folder if available (`State::load`), otherwise default state is created.
4. CLI settings/runtime defaults are loaded and synchronized with tracing filters.

### 2. Process Generation Flow
1. `generate` command builds `ProcessDefinition` (from syntax or graph import).
2. `State::generate_integrand(s)` creates a generation thread pool.
3. `ProcessList::preprocess` delegates to amplitude/cross-section preprocessors.
4. `ProcessList::generate_integrands` builds `ProcessIntegrand` instances from preprocessed graphs.
5. Optional compile/export steps persist compiled evaluator artifacts and DOT/standalone outputs.
6. Each generated integrand now embeds its frozen f64 backend choice in `integrand.bin`:
   - `eager`
   - `symjit`
   - external `c++` / `assembly` plus the external compile options used
7. Runtime-only evaluator backends are then activated from that frozen metadata:
   - eager uses the saved eager evaluator directly
   - symjit is rebuilt after generation/load from the saved Symbolica evaluator
   - external compiled backends load their saved shared-library artifacts lazily
   - if external loading fails and startup globals explicitly opt into symjit,
     GammaLoop falls back to symjit for that integrand and logs it

### 3. Evaluation and Integration Flow
1. Commands (`inspect`, `evaluate`, `integrate`) resolve process + integrand references.
2. Integrand is warmed up (`ProcessIntegrand::warm_up`) to initialize rotations and caches.
   Every evaluator receives the auxiliary numerator sampling scale `M` as an
   input, defaulting to one. Warm-up rejects a zero runtime scale for every
   amplitude or cross-section integrand, including evaluators that do not use
   `M`. No source-expression scan or serialized usage flag is needed. This
   enforces the EMR-only `a*M` contract above; `M` is not an LMB coordinate.
3. Sampling path parameterizes points and evaluates graph terms.
4. Stability checks may escalate precision (`f64 -> f128 -> arbitrary`) and rotate kinematics.
5. Process graph evaluation returns a rich `GraphEvaluationResult<T>` rather
   than only a complex weight. This carries:
   - the graph contribution
   - grouped generated events
   - event-processing timing
   - generated / accepted event counts
6. A shared process-layer prepared-event helper is used by both cross-section
   and amplitude graph evaluation:
   - it decides whether an event must be built for the current rotation
   - runs selectors immediately on that candidate event
   - retains only accepted identity-rotation events when buffering is needed
   - reports generated / accepted event counts and event-processing timing
   Failed selector events zero the local contribution and do not survive into
   the final retained result.
7. Stability selection still compares only the complex graph weight, but the
   retained branch also carries the final grouped event payload.
8. The final `EvaluationResult` contains:
   - the stable `integrand_result`, before any parameterization Jacobian is
     applied
   - the top-level `parameterization_jacobian` when the sample came from
     x-space parameterization (`None` for direct momentum-space evaluation)
   - the separate `integrator_weight`, i.e. the Monte Carlo/grid weight only
   - grouped accepted events
   - event-processing timing in `evaluation_metadata`
   - generated / accepted event counts in `evaluation_metadata`
   - ordered per-level stability results, each with relative-accuracy and total
     time spent in that stability level
   - evaluation metadata
9. Observable filling happens only from the final stable `EvaluationResult`
   payload. Unstable branches do not contribute events or observables.
10. `havana_integrate` runs iterative Monte Carlo updates, merges worker-local
    observable accumulators, and writes integration artifacts.

### 3.1 Differential event-processing runtime

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

The runtime-only cache is intentionally not serialized in binary state dumps.
It is wrapped in a no-op `Encode` / default-empty `Decode` container and is
always rebuilt from settings at warm-up time.

The same pattern is now also used for evaluator execution backends:

- per-integrand frozen backend metadata lives inside `integrand.bin`
- loaded external `.so` evaluators are runtime-only
- symjit evaluators are runtime-only
- the saved portable representation remains centered on eager Symbolica
  evaluators

### 3.2 Differential event model

Events are stored in precision-homogeneous containers:

- `GenericEvent<T>`
- `GenericEventGroup<T>`
- `GenericEventGroupList<T>`

This avoids enum-based mixed-precision layouts that would otherwise force
storage to be sized by the largest precision variant.

Each event carries:

- kinematics
- cut metadata
- one fully normalized event `weight`
- `additional_weights`

Observable-specific entry reweighting remains internal to the observable
runtime; it is not stored on the event.

`additional_weights` is a generic `BTreeMap` keyed by lightweight identifiers
such as:

- `FullMultiplicativeFactor`
- `Original`
- `ThresholdCounterterm { subset_index }`

When threshold-counterterm weights are stored, they are stored with the same sign
with which they contribute to the final event weight. This means the fully
normalized event weight is reconstructed as:

`(Original + sum(ThresholdCounterterm { subset_index })) * FullMultiplicativeFactor`.

This is populated only when
`settings.general.store_additional_weights_in_event = true`.

### 3.3 Event generation policy

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
  - selector-only evaluations still build temporary events for selector checks
  - selector+observable evaluations still build temporary events and buffer
    accepted identity-rotation events long enough to fill observables
  - those temporary events are cleared from returned sample results and
    integration event outputs afterwards
- `general.generate_events = true`
  - events are built, buffered, and returned even when no selector or
    observable is configured

The standard Rust/Python sample-evaluation helpers may temporarily override the
returned-event behavior per call, but internal event generation still follows
this same policy after that override is applied.

Accepted LU events are retained as grouped event lists:

- one `EventGroup` per graph-group evaluation
- all accepted cuts from graphs sharing the same `group_id` are merged into the
  same retained event group
- each event stores its concrete `graph_id` and `cut_id`
- incoming PDGs are sourced from the initial-state cut edges

This preserves the correlation structure needed by downstream observables.

### 3.4 Observable architecture

Observable definitions are shared between selector and histogram use-cases:

- one quantity definition extracts zero or more `ObservableEntry<T>` values from
  an event
- selectors consume those entries on a single event
- histogram observables consume grouped events and fill per-group contributions

Histogram accumulation is based on:

- a histogram-level `sample_count`
- sparse per-bin sufficient statistics:
  - `entry_count`
  - `sum_weights`
  - `sum_weights_squared`
  - `mitigated_fill_count`
- `HistogramAccumulatorState` / `ObservableAccumulatorBundle` for mergeable
  runtime state
- `ObservableSnapshotBundle` for JSON / HwU output and API responses

Each Monte Carlo sample contributes exactly one statistical sample to a
histogram. Contributions from all retained event groups are summed first per
bin, and untouched bins are handled implicitly through the histogram-level
`sample_count` rather than explicit zero fills.

Histogram snapshots therefore remain fully mergeable and re-constructible into
live accumulator state. The Rust-facing histogram API exposes `merge(...)`,
`merge_in_place(...)`, and `rebin(...)`.

Current built-in quantities include:

- `particle` with `computation = scalar | count | pair`
- `jet` with `computation = scalar | count | pair`
- `afb`
- `cross_section`

Current pair quantities include `DeltaR`. Current scalar projections include
`E`, `CosTheta`, `PT`, `y`, `eta`, `Px`, `Py`, `Pz`, and `Mass`.

Object and pair quantities are ordered before `entry_selection` is applied.
The ordering contract is configurable per quantity through:

- `ordering = PT | Energy | AbsRapidity | Quantity`
- `order = Ascending | Descending`

Family defaults are:

- particle scalar quantities: `ordering = Quantity`
- jet scalar quantities: `ordering = PT`
- pair quantities: `ordering = Quantity`

`leading_only` and `nth_only` therefore operate on this explicit ordered list,
not on incidental source/event enumeration order.

Jet clustering is implemented natively in Rust with IRC-safe algorithms:

- `kt`
- `cambridge_aachen`
- `anti_kt`

and is validated against `fjcore` in tests.

## Persistence and Runtime Artifacts
Primary persisted state lives under `gammaloop_state/` (default):
- `state_manifest.toml` (state schema/version marker)
- `model.json`, `model_parameters.json`
- `symbolica_state.bin`
- `processes/` (amplitudes/cross_sections + integrands)
- `run.toml`
- `default_runtime_settings.toml`
- `global_settings.toml`
- `logs/`

For local experimentation, prefer an isolated path such as `.local/scratch/<run>/gammaloop_state` to keep repository root output minimal.

The persistence model is file-system based and intentionally human-editable for settings/run cards, mixed with binary artifacts for performance-heavy data.

### Persistence Compatibility Contract
- State format is versioned with `state_manifest.toml` (`version = 6` currently).
- Version 6 removes obsolete deferred-integrand fields from saved amplitudes
  and cut integrands. Version 5 added component-local generated-CFF ownership
  and prefactor metadata; version 4 added the typed CFF core global-prefactor
  sign to the version-3 layout. These changes affect positional bincode data,
  so version-3, version-4 and version-5 states must be regenerated rather than
  relabeled.
- State loading and direct overwrite both require exactly the current manifest version; older states must be regenerated, and states from newer binaries require a newer GammaLoop binary.
- A missing manifest denotes an unmanifested folder rather than a legacy state and is never loaded as saved state.
- Process settings history now uses `settings_history.toml` consistently; loader still accepts legacy `settings_history.yaml` for backward compatibility and migration.

## Configuration Architecture
Configuration is split into:
- Global settings (`GlobalSettings`): generation, logging directives, parallelism.
- Runtime settings (`RuntimeSettings`): kinematics, integrator behavior, sampling, stability, subtraction, and differential event-processing configuration.

Differential runtime configuration currently includes:

- `general.generate_events`
- `general.store_additional_weights_in_event`
- named `quantities`
- named `selectors`
- named `observables`
- `integrator.observables_output`

This split is good: generation-time and evaluation-time concerns are separated and independently serializable.

## Concurrency Model
Concurrency is explicit and use-case scoped:
- Generation and compile thread pools use configurable thread counts.
- Integrator parallelism is controlled via runtime/global settings.
- Some loops over processes/integrands remain sequential at orchestration level while heavy operations inside are parallelized.

For the differential pipeline:

- each worker/integrand clone owns its own mutable observable accumulator state
- local worker results are merged back into the master integrand after chunk
  evaluation
- distributed batch mode can return either grouped events or serializable
  observable accumulator bundles
- observable snapshots are written only from the merged master state

## Testing Architecture
Testing is organized into layers:
- Unit tests in modules across `crates/gammalooprs/src/`.
- Integration test crate: `tests/` with reusable harness (`tests/src/lib.rs`).
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
  - real generated-process end-to-end coverage for Rust `evaluate_sample(s)`
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



## Architectural Strengths
- Clear separation between application shell (`gammaloop-api`) and domain core (`gammalooprs`).
- Rich, serializable settings model with schema generation hooks.
- Strong process abstraction (`Process`, `ProcessCollection`, `ProcessIntegrand`) that supports amplitude and cross-section workflows.
- Robust stability pipeline with multi-precision escalation and rotation checks.
- Differential event-processing is now integrated into the core evaluation path
  without compromising the stable-weight selection logic.
- Observable accumulation is mergeable across rayon workers and distributed
  batch outputs.

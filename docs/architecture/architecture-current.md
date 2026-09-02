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

All CFF power and capacity questions, including numerator, repeated-channel,
and finite-pole powers, are expressed solely in physical EMR/source-edge
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
up to sign to validate denominator equality and expose repeated CFF channels;
they never infer endpoints or merge physical owners. The only nearby rank solve
selects a unique `+/-` routing sign modulo the opposite source domain on those
already fixed endpoints. There is no incidence/Kirchhoff reconstruction and no
external-balance synthesis. The raw `+/-Q` sign remains available to the
numerator mapper, while a post-construction Graphica pass canonically relabels
nodes and exact edges for deterministic cache keys.

Projection first plans every genuinely outer additive Taylor term, then
generates once per canonical topology. Every term keeps its own factorized
minimax assignment. A shared generator envelope takes the maximum **total**
degree within each repeated algebraic energy channel and redistributes that
total by the same minimax rule; non-repeated bounds use componentwise maxima.
The real degree-one triangle regression retains its collected coefficient as
one common-denominator source with UV-owner multiplicities `(2,1,2)`. Positive
typed denominator factors stay in its factorized numerator, and generalized
CFF produces the lower sectors internally. It therefore uses one canonical
generator entry and agrees with an explicitly reduced test-only oracle.
Non-vacuum exact sources retain pure-external boundaries as explicit
source-crown hedges. Future on-shell two-point insertions such as `(m,0,0,0)`
require an explicit fixed-boundary payload, not topology reconstruction.
This exact-source reconstruction and minimax-dispatch machinery is exclusive to
the projected local-4D route.
The full design, sign argument for `D(Q)=D(-Q)`, and concrete production fixtures
are documented in
[`exact-powered-denominator-cff-lifting.md`](exact-powered-denominator-cff-lifting.md).

The production numerator remains factorized. Degree analysis traverses its
factors without expanding them, each UV step attaches only newly owned factors,
and final assembly attaches outside and global factors exactly once. For higher
power projection, the term parser splits only the completed expression's outer
Taylor sum when its addends carry separate denominator topologies. Nested
numerator sums and positive typed `GS.den` factors remain factorized and are
never cancelled against denominator occurrences upstream; generalized CFF owns
the resulting pinches and lower sectors. For higher powers,
interpolation may replace an EMR energy by `a*M`, where `a` is a signed
integer and `M` is the common auxiliary CFF numerator-sampling scale. This is an
EMR substitution, never an LMB rewrite. Only finalized evaluators that use `M`
require `M != 0`, and the physical result is invariant under changing its
nonzero value.

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
The empty UV forest is the ordinary factorized production root in both local-UV
routes; the expanded-4D setting changes only proper, nonempty UV nodes.
The shared CFF core also returns its connected-loop and pure
duplicate-denominator global sign as typed metadata. GammaLoop consumes that
bridge exactly once for root, reduced, and exact production CFF sources,
cancelling the shared-core-local uniform convention and retaining GammaLoop's
established complete-integrand convention.
The NLO acceptance layer exercises these production boundaries in
orientation-local direct 3D, explicit-sum direct 3D, and projected local 4D.
The retained 2026-08-31 10x campaign completes all four physical DD/TT
acceptances (`4/4`). Pulls below are signed differences from the published
target in units of the reported Monte Carlo error; ratio pulls include the LO
uncertainty.

| Acceptance | 10x LO result | 10x NLO result | Graph and ratio evidence |
| --- | --- | --- | --- |
| direct `gamma* -> d d~` | `0.5068703962 +/- 0.0025987972` (`+1.449 sigma`) | `0.01966009810 +/- 0.00053595339` (`+1.424 sigma`) | `GL0=-0.03132123586 +/- 0.00023922726` (`+0.729 sigma`), `GL2=+0.05112479005 +/- 0.00046213299` (`+1.584 sigma`); `alpha_s/pi` pull `+1.141 sigma` |
| converted `e+e- -> gamma* -> d d~` | `0.1950499744 +/- 0.0010326753 pb` (`+1.479 sigma`) | `0.007824189766 +/- 0.000340601513 pb` (`+1.630 sigma`) | signed MC components `GL0=-0.01996339254 +/- 0.00015479207`, `GL2=+0.02786745983 +/- 0.00028425324`; they have no separate published targets; `alpha_s/pi` pull `+1.453 sigma` |
| direct `gamma* -> t t~` | `2.901968994 +/- 0.015639978` (`+1.641 sigma`) | `0.2079169992 +/- 0.0042953541` (`+1.489 sigma`) | `GL0=-0.1443600613 +/- 0.0035809931` (`+0.669 sigma`), `GL2=+0.3522770605 +/- 0.0023720361` (`+1.687 sigma`); paper-ratio pull `+1.037 sigma` |
| converted `e+e- -> gamma* -> t t~` | `0.3307052414 +/- 0.0018004843 pb` (`+1.603 sigma`) | `0.02356890542 +/- 0.00056839205 pb` (`+1.058 sigma`) | summed-graph integration has no persisted GL0/GL2 rows; paper-ratio pull `+0.685 sigma` |

All LO integrations used 100,000 samples. The direct DD NLO used 400,000
samples, the direct TT graph rows used 400,000 each, and each converted NLO
central slot used 200,000. The converted DD acceptance also retains its
pointwise route, scale-law, and EMR-bound checks; its physical graph components
are closure diagnostics rather than separately published observables.

The scalar local-equivalence matrix is generated from the scalar model rather
than from hand-built graph data. Its unit-numerator lanes remain unchanged after
generation, companion probes use only Feynman-rule-local edge factors, and
there is no graph-specific production branch. The scalar LU 15-case matrix
enables local UV, integrated UV, and threshold counterterms while comparing
orientation-local 3D, explicit-sum 3D, and projected local 4D, including native
Arb checks. As retained pre-reversal evidence, its 2026-08-31 `dev-optim` /
`test_gammaloop` rerun passed `15/15`, with `235` skipped, in `70.438 s`. The
matrix must be rerun after restoring the post-CFF direct-local3D construction
before it is a current merge-readiness gate. For four near-zero cases, the
authorized f64-input comparison uses the `1e-14` unit-scale fallback; the
Arb-to-Arb comparison still runs and reports non-scaling. This is test-oracle
handling only and required no production change. The curated suite selects a
focused DOD0/1/2
orientation-local bubble regression and all six base scalar-matrix graphs
(GL00, GL02, GL04, GL08, GL09, and GL24); each base graph now also invokes
per-orientation profiling for the localized direct-local3D route. Their current
post-restoration rerun is pending.

UV profiling defaults to `only-divergent`: every expected cycle union with
DOD >= 0 is tested using the generation LMB when suitable, otherwise the first
suitable basis in the deterministically sorted complete LMB list. The
exhaustive `all` mode is opt-in. Amplitude and LU inputs share this behavior,
with graph and Cutkosky-cut selectors for LU profiling and a colored final
failure summary. Current end-to-end CLI coverage is `2/2` and lower-level
API/unit coverage is `15/15`. Per-orientation profiling is defined only for the
orientation-parametric, localized direct-local3D mode. Selector-free explicit-
sum direct local3D and projected local4D are summed representations and reject
that request. The shared CFF crate is currently `98/100`; its
only failures are the powered-pole fixture contract and the inherited
mixed-theta origin expectation. Both await test-only cleanup and are not a
demonstrated production-value mismatch.

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
   Each finalized evaluator persists whether its source expressions (including
   deferred function bodies) contain the auxiliary numerator sampling scale
   `M`; warm-up rejects a zero runtime scale only when at least one evaluator
   in the selected amplitude or cross-section integrand actually uses it. This
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
- State format is versioned with `state_manifest.toml` (`version = 4` currently).
- Version 4 persists the typed CFF core global-prefactor sign. Because this
  changes the positional bincode layout of generated three-dimensional
  expressions, version-3 states must be regenerated rather than relabeled.
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

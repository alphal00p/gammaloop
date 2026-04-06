# Discrete Histograms

## Status

- Date: 2026-04-04
- Phase: implementation completed and verified at compile/test level
- Implementation status:
  - tagged continuous/discrete histogram settings and snapshots implemented
  - discrete metadata quantities and event metadata propagation implemented
  - selector `active` / `discrete_range` and observable named selections implemented
  - Rust and Python histogram accumulator APIs exposed
  - CLI/settings templates and schema/completion support updated
- checked-in example cards for `cross_section` / `jet_count` migrated to discrete definitions
- `examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml` now includes advanced discrete `graph_id`, per-graph `orientation_id`, and per-graph `lmb_channel_id` observables driven by inactive named graph selectors and Monte Carlo sampling over graphs/orientations/LMB channels
  - verification complete with `cargo fmt`, `cargo check`, `cargo check --all-targets`, focused Rust/API tests, and `cargo clippy --all-targets`
  - remaining clippy warnings are pre-existing repository baseline outside this feature slice

## Delivered Implementation

- `CutInfo` now keeps one shared numeric metadata model for runtime and output with:
  - `graph_id`
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
  - `lmb_channel_edge_ids`
- Added discrete observable quantities:
  - `graph_id`
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
- Added selector/runtime changes:
  - selector `active: bool`
  - inclusive `discrete_range`
  - observable `selections: Vec<String>` as named selector references
  - active selectors gate runtime event retention, while observable selections gate only filling
- Histogram infrastructure now supports continuous and discrete variants with:
  - stable `bin_id`
  - optional labels
  - ordering strategies
  - merge by `bin_id`
  - no-op discrete `rebin()`
  - HwU discrete serialization as `[k, k + 1)`
- Auto discrete domain / label resolution is implemented for:
  - graph ids
  - graph-group ids
  - orientation ids
  - LMB channel ids
  using singleton graph/group selector context or single-group integrands
- Count-like quantities and `cross_section` now fill through the discrete path.
- Standalone histogram accumulation is exposed through pyo3 with fill, merge, snapshot, rescale, rebin, and ordering APIs.
- Observable / selector serde and schema stay flat at the card level, while the internal histogram model is tagged continuous/discrete.

## Verification

- `cargo fmt`
- `cargo check`
- `cargo check --all-targets`
- `cargo test -p gammalooprs observable_settings_infer_discrete_histogram_from_domain -- --nocapture`
- `cargo test -p gammalooprs selector_settings_serialize_flat_discrete_range -- --nocapture`
- `cargo test -p gammalooprs selector_settings_deserialize_flat_discrete_range -- --nocapture`
- `cargo test -p gammalooprs observable_settings_schema_is_flat -- --nocapture`
- `cargo test -p gammalooprs selector_settings_schema_is_flat -- --nocapture`
- `cargo test -p gammalooprs runtime_test_serialize_deserialize -- --nocapture`
- `cargo test -p gammalooprs discrete_histogram_ordering_changes_presentation_only -- --nocapture`
- `cargo test -p gammaloop-api parse_set_process_update_observable -- --nocapture`
- `cargo test -p gammaloop-api process_observable_accepts_histogram_title_and_type_description -- --nocapture`
- `cargo clippy --all-targets`
- `./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate -c "inspect -p epem_a_tth -i LO -x 0.1 0.2 0.3 0.4 0.5 0.6; quit -o"`
- `./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate set_model_parameters set_differential_quantities set_differential_selectors set_differential_observables -c "set process -p epem_a_tth -i LO string ' [sampling] graphs = \"monte_carlo\" orientations = \"monte_carlo\" lmb_multichanneling = true lmb_channels = \"monte_carlo\" '; inspect -p epem_a_tth -i LO -d 0 0 0 -x 0.1 0.2 0.3 0.4 0.5 0.6; quit -o"`
- short LO integration cross-checks at `5k` and `20k` samples for:
  - summed sampling
  - Monte Carlo graph/orientation/LMB sampling
  - Monte Carlo graph/orientation/LMB sampling plus the new observables/selectors

## Notes

- The implementation keeps one shared numeric event-metadata model. No public/internal split was introduced.
- Returned histogram bundles use concrete labels only; auto-label strategies resolve during runtime initialization, not during hot filling.
- Current repository-wide clippy warnings still exist outside this feature area; the discrete-histogram work was reduced to zero new local clippy warnings.
- The LO `epem_a_ttxh` example card now demonstrates the intended advanced configuration pattern:
  - `graph_id` histogram with `GraphIds` + `GraphName`
  - per-graph `orientation_id` histograms using named inactive `graph_id` selectors plus `OrientationIds` + `Orientation`
  - per-graph `lmb_channel_id` histograms using named inactive `graph_id` selectors plus `LmbChannelIds` + `LmbChannelEdgeIds`
  - Monte Carlo sampling over graphs, orientations, and LMB channels
- Investigation follow-up for the LO `epem_a_ttxh` example:
  - the earlier `inspect` regression was due to `inspect` not forcing `generate_events`; this is fixed
  - the new inactive named selectors behave as intended:
    - they do not participate in runtime event rejection
    - they do gate only the corresponding observables
  - with discrete Monte Carlo sampling active, `inspect -d 0 0 0 ...` shows `graph_id = 0`, `orientation_id = 0`, and `lmb_channel_id = 0` in event metadata, and only the `graph_0_*` discrete histograms are filled
  - integration with and without the new observables/selectors gives the same estimate for the same seed and sample count
  - compared to the old summed sampling, Monte Carlo sampling over graph/orientation/LMB choices is unbiased in the tested runs but has visibly larger variance, so convergence is slower/noisier for this example

## Approved Implementation Scope

- Keep one shared numeric `CutInfo` for both runtime and returned output. No internal/output split.
- Add optional metadata fields:
  - `graph_group_id: Option<usize>`
  - `orientation_id: Option<usize>`
  - `lmb_channel_id: Option<usize>`
  while keeping `graph_id` and `lmb_channel_edge_ids`.
- Expose four discrete metadata-backed quantities:
  - `graph_id`
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
- Replace the flat histogram definition with a tagged continuous/discrete model.
- Discrete domains:
  - `ExplicitRange { min, max }`
  - `SingleBin`
  - `GraphIds`
  - `GraphGroupIds`
  - `OrientationIds`
  - `LmbChannelIds`
- Discrete labels:
  - `Custom`
  - `BinId`
  - `GraphName`
  - `GraphGroupMasterName`
  - `Orientation`
  - `LmbChannelEdgeIds`
- Discrete ordering:
  - `AscendingBinId` default
  - `ValueDescending`
  - `AbsValueDescending`
- Add `change_bin_ordering(...)` in Rust and Python APIs. Ordering is presentation-only; stable `bin_id` drives fill and merge semantics.
- Add selector `active: bool` and a new inclusive `DiscreteRange { min, max }` selector mode.
- Add `ObservableSettings.selections: Vec<String>` as named selector references only.
- Runtime event filtering uses only active selectors. Observable-local selections can reference both active and inactive selectors and gate only observable filling.
- Group-local auto domain and auto-label resolution for orientations and LMB channels may use:
  - a singleton referenced selector on `graph_group_id`
  - a singleton referenced selector on `graph_id` mapped to its graph group
  - an integrand with exactly one graph group
- `graph_group_id` must always be populated whenever `graph_id` is known.
- `orientation_id` and `lmb_channel_id` are populated only when the normalized sample fixes them.
- Migrate count-like observables and `cross_section` to the discrete path. `cross_section` uses `SingleBin` with default label `total cross-section`.
- Keep discrete HwU output as `[k, k + 1)`, disable misbinning for discrete histograms, and make discrete `rebin()` / contiguous-bin merge no-ops.
- Expose the standalone histogram core through pyo3, and update CLI completions/schema/templates for the new settings.

## Objectives

- Add a first-class discrete histogram variant alongside the current continuous histogram path.
- Keep the observable-facing API close to the current histogram API.
- Preserve current HwU output conventions for discrete-like observables.
- Expose the histogram infrastructure itself through pyo3, not only immutable evaluation snapshots.
- Add discrete metadata-driven quantities such as `graph_group_id`, `orientation_id`, and `lmb_channel_id`.
- Extend event metadata with optional discrete IDs so metadata-driven quantities can skip filling cleanly when unavailable.
- Keep returned event metadata numeric-only and shared with the runtime path.
- Allow observable-local selections so histogram filling can be constrained without relying only on global selectors.
- Migrate existing jet-count and cross-section quantities plus existing TOML cards to the new discrete model.
- Preserve and update CLI/settings completion for the new quantity and histogram types.

## Current Findings

- The current histogram stack is continuous-only end to end:
  - `HistogramSettings`
  - `HistogramAccumulatorState`
  - `HistogramSnapshot`
  - `HistogramObservable`
  - `ObservableSnapshotBundle` / `ObservableAccumulatorBundle`
- Count-like observables are only special at the quantity layer. `ResolvedQuantityComputation::Count` disables misbinning mitigation, but still emits a scalar value that is binned through the continuous histogram path.
- The low-level filling logic is not currently a standalone histogram API. It is embedded inside `HistogramObservable::build_group_contributions`, `mitigate_group_misbinning`, and `flush_sample_contributions`.
- The Python API currently exposes observable snapshots in evaluation results, but not standalone histogram creation/fill/merge/rebin/write APIs.
- The Python API also currently exposes:
  - `PyCutInfo.graph_id: usize`
  - a continuous-only `PyHistogramSnapshot`
  so both event metadata and histogram snapshots will need a public shape change.
- The current selector compiler is already reusable:
  - global selectors compile into quantity-based event predicates
  - observables are compiled separately
  so observable-local selections can likely reuse the same compiled predicate machinery instead of inventing a second quantity-evaluation path.
- Settings completion is schema/serialization-driven through:
  - `serialize_settings_with_defaults(...)`
  - `serialize_schema(...)`
  - `SHOWDEFAULTS`
  so changes to runtime settings types and their `JsonSchema`/serde shape will directly affect completion.
- Event metadata already carries:
  - `graph_id` on `CutInfo`
  - `lmb_channel_edge_ids` on `CutInfo`
- Event metadata does not currently carry:
  - `graph_group_id`
  - `orientation_id`
  - a compact `lmb_channel_id`
- Metadata origin is already available in the evaluation flow:
  - `graph_id` is assigned centrally in `evaluate_graph_term(...)`
  - graph-group id is available at group-evaluation level
  - orientation is stored on `MomentumSample.sample.orientation`
  - channel id is available during event construction, but only edge ids are currently persisted
- `evaluate_sample(s)` already accepts graph/group/orientation/channel information through more than one route:
  - x-space `discrete_dims`
  - momentum-space `discrete_dims`
  - explicit momentum-space `graph_name`
  - explicit momentum-space `orientation`
  - explicit momentum-space `channel_id` once a graph group is selected
  These routes converge through `GammaLoopSample` / `DiscreteGraphSample`, so metadata attachment should happen after that normalization step rather than in each frontend entrypoint.
- Existing display helpers already define the user-facing formatting conventions we should reuse:
  - orientations: uncolored `+`, `-`, `0` signature strings
  - LMB edge lists: `(e1,e2,...)`
  - graph groups: master-graph naming already exists in integrand/group helpers
- Important semantic finding:
  - `graph_group_id` is global
  - `orientation_id` is local to a graph group
  - `lmb_channel_id` is local to a graph group
  so the automatic label strategies `Orientation` and `LmbChannelEdgeIds` are not globally meaningful unless the histogram is effectively tied to one graph group, or the mapping is proven identical across all contributing graph groups.
- Important output-model finding:
  - `EvaluationResultOutput` currently serializes `EventGroupList` directly
  - `EventGroupList` uses the same `CutInfo` type as the internal event-processing path
  so a true public/internal metadata split is not just a Python change; it requires dedicated output event structs and conversion in `EvaluationResult::into_output(...)` and the precise equivalents.

## Confirmed Decisions

- Clean break:
  - existing jet count and cross-section move to the new discrete model
  - existing TOML cards found in the repo should be migrated
- JSON should use a tagged continuous/discrete shape.
- Discrete histogram `rebin(...)` should leave the histogram unchanged.
- HwU discrete convention is locked to numeric bins `[k, k + 1)`.
- Missing optional metadata means:
  - event metadata field stays `None`
  - quantity extraction emits no entry and therefore no fill
- `lmb_channel_id` means exactly the discrete sampling value used for LMB multichanneling.
- What was previously called `graph_id` for the new discrete quantity should actually be `graph_group_id`, matching the discrete graph-sampling dimension semantics.
- Returned event metadata should expose `graph_name` rather than a numeric `graph_id`.
- Returned histogram bundles from `evaluate_sample(s)` should only contain concrete labels, not unresolved auto-label strategy enums.
- Add a `BinID` automatic label strategy producing `#<bin_id>`.
- If `orientation`, `graph_group`, or `lmb_channel_id` is supplied to `evaluate_sample(s)` through any supported route, it must end up both in event metadata and in the corresponding discrete histogram quantities.
- Default discrete bin ordering should be stable ascending bin id.
- The API should allow changing discrete bin ordering later, with fill/merge compatibility always based on stable `bin_id`, not rendered order.
- The evaluation path is hot code:
  - no per-fill label resolution
  - no per-fill string allocation
  - attach compact metadata once and reuse it everywhere

## Proposed Direction

### 1. Split histogram core from observable-specific extraction

- Extract a standalone histogram engine that owns:
  - histogram definition/layout
  - grouped sample filling
  - merge/update/snapshot/restore
  - JSON and HwU writing
- Keep `HistogramObservable` as a thin adapter:
  - build observable entries from events
  - hand grouped entries to the standalone histogram engine
- This is required for clean pyo3 exposure. Otherwise the Python API would have to go through the whole observable runtime just to fill a histogram.

### 2. Introduce explicit histogram layout variants

- Replace the implicit continuous-only layout with an explicit axis/layout enum, conceptually:
  - `Continuous { x_min, x_max, n_bins, log_x_axis }`
  - `Discrete { min, max, labels, sort }`
- Keep shared histogram-level fields common:
  - `title`
  - `type_description`
  - `phase`
  - `value_transform`
  - `sample_count`
  - `log_y_axis`
  - shared per-bin statistics
- Discrete layout rules:
  - `min` and `max` define an inclusive integer range of in-range bin ids
  - each in-range bin covers exactly one integer unit
  - underflow and overflow bins remain explicit
  - misbinning mitigation is always unsupported
  - `bin_id` stays separate from display order
  - bin labels are optional metadata attached to bin ids, not bin positions
  - discrete `rebin(...)` is a no-op returning the histogram unchanged

### 2b. Make discrete bin ordering explicit

- Discrete histograms need two separate concepts:
  - stable `bin_id` used for filling, merging, compatibility checks, and label lookup
  - current ordered bin list used for JSON/API rendering
- Recommended direction:
  - keep internal per-bin storage keyed by `bin_id`
  - maintain a separate ordered list or sortable view for output
- Add configurable sort strategies in histogram definition and API, for example:
  - `AscendingBinId` as the default
  - `ValueDescending`
  - `AbsValueDescending`
- Sorting must be changeable after construction through the Rust and Python APIs without changing histogram identity or merge semantics.
- Expose a dedicated API such as `change_bin_ordering(...)` for this purpose.

### 3. Introduce explicit observable coordinates

- Refactor `ObservableEntry<T>` so its coordinate is no longer always continuous.
- Recommended direction:
  - `ObservableCoordinate<T> = Continuous(F<T>) | Discrete(isize)`
  - `ObservableEntry<T> { coordinate, weight_modifier }`
- Continuous quantities continue to emit `Continuous(...)`.
- Discrete quantities emit `Discrete(...)`.
- This gives a real `isize` path internally instead of pretending discrete values are floats.

### 4. Add discrete-native quantities

- Add hardcoded quantity variants for:
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
- Migrate naturally discrete existing quantities to the new coordinate model:
  - `cross_section`
  - count-based particle/jet quantities
- Metadata-driven quantities should return zero entries when metadata is unavailable.
  - This is cleaner than inventing sentinel values.
  - It naturally means the histogram is not filled for that event/group.

### 5. Normalize discrete evaluation context once

- Introduce a small internal context extracted from the normalized evaluation sample, conceptually:
  - `graph_group_id: Option<usize>`
  - `orientation_id: Option<usize>`
  - `lmb_channel_id: Option<usize>`
- Build this context from `GammaLoopSample` / `DiscreteGraphSample` rather than from raw API arguments.
- This guarantees that all routes into evaluation behave the same:
  - x-space `discrete_dims`
  - momentum-space `discrete_dims`
  - explicit `graph_name`
  - explicit `orientation`
  - any future path that also lowers into the same sample types
- Populate metadata from the evaluation context once and then let quantities read it without extra branching or reconstruction.

### 5b. Split internal metadata from output-facing metadata

- Recommended direction:
  - keep the hot internal event metadata compact and numeric
  - enrich it to output-facing `graph_name` only when building returned/public event results
- Internally, the event path should still carry:
  - `graph_id`
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
  - `lmb_channel_edge_ids`
- Public evaluation output should expose:
  - `graph_name`
  - `graph_group_id`
  - `orientation_id`
  - `lmb_channel_id`
  - `lmb_channel_edge_ids`
- This avoids per-event owned-string allocation inside the hot generation/filling path while still giving the caller the more useful public metadata.
- Complexity assessment:
  - this is moderate, not tiny
  - it requires dedicated output event/cut-info structs and conversion in the evaluation-output layer
  - it does not require changing histogram filling logic itself
- Current recommendation:
  - keep the split
  - the scope is still contained enough to be worth it if we want to avoid repeated graph-name allocation in the hot event path

### 5c. Add observable-local selections and use them to resolve discrete context

- Add a new observable field:
  - `selections: Vec<ObservableSelectionSettings>`
  - default `[]`
- Recommended observable-selection variants:
  - `Selector(SelectorSettings)` using an inline selector definition
  - `GraphGroupIDEqual(usize)`
  - `OrientationIDEqual(usize)`
  - `LMBChannelIDEqual(usize)`
- Semantics:
  - all selections are ANDed
  - they are evaluated per event before the observable extracts entries
  - global selectors still run first and determine whether the event survives at all
- Implementation direction:
  - reuse the compiled selector predicate machinery for the inline-selector variant
  - implement the three metadata-equality variants as direct integer checks with no quantity evaluation
- These observable-local selections then become the source of context for auto-resolved discrete labels and auto-resolved discrete bin ranges.

### 5d. Add automatic discrete bin-label strategies, but resolve them once

- The observable/settings layer should allow:
  - `Custom(["label1", "label2", ...])`
  - `BinID`
  - `Orientation`
  - `GraphsGroupMasterName`
  - `LmbChannelEdgeIds`
- Strategy semantics:
  - `Custom`: positional labels for the inclusive integer interval `[min:max]`
  - `BinID`: generate `#<bin_id>`
  - `Orientation`: derive the uncolored `+-0` signature string for the corresponding orientation id
  - `GraphsGroupMasterName`: derive the master graph name for the corresponding graph-group id
  - `LmbChannelEdgeIds`: derive the sorted edge-id list string for the corresponding LMB channel id
- Process-aware label strategies must be resolved exactly once:
  - preferably during observable compilation/runtime initialization when the integrand is known
  - or, if absolutely necessary, when building the final snapshot bundle
- The snapshot and Python-facing histogram types should only ever contain concrete labels.
- Standalone pyo3 histogram construction should therefore expose only concrete labels; `BinID` can be resolved immediately into concrete labels too.

### 5e. Add auto-resolved discrete bin ranges

- Instead of only a literal inclusive integer range, the discrete histogram definition should support a range enum, conceptually:
  - `Explicit { min, max }`
  - `GraphGroupIDs`
  - `OrientationIDs`
  - `LMBChannelIDs`
- Resolution rules:
  - `GraphGroupIDs` resolves globally to `[0, graph_group_count - 1]`
  - `OrientationIDs` resolves to `[0, n_orientations(group) - 1]`
  - `LMBChannelIDs` resolves to `[0, n_channels(group) - 1]`
- Because discrete ranges are inclusive, any count-based auto range resolves to `count - 1`, not `count`.
- Resolution context for `OrientationIDs` and `LMBChannelIDs`:
  - use `GraphGroupIDEqual(...)` from the observable selections when present
  - otherwise, if the integrand has exactly one graph group, use that implicitly
  - otherwise fail early with a clear error
- The same resolution context should be used for `Orientation` and `LmbChannelEdgeIds` auto labels.

### 6. Keep HwU output identical in practice

- The discrete histogram variant should still write numeric bin edges to HwU.
- Recommended convention for exact compatibility with current count-style histograms:
  - integer bin `k` is written as `[k, k + 1)`
  - a one-bin cross-section histogram remains `[0, 1)`
- Bin labels remain JSON/Python metadata only and do not change HwU format.

### 7. Expose the histogram stack through pyo3

- Export standalone histogram types, not only evaluation-result snapshots.
- Target Python surface:
  - histogram definition/layout objects
  - accumulator/state objects
  - snapshot objects
  - bundle objects
  - methods for:
    - construction
    - grouped filling
    - merge
    - update/finalize
    - snapshot/restore
    - rebin / merge contiguous bins
    - change discrete bin sort strategy
    - JSON write/read
    - HwU write
- Keep evaluation-result observable snapshots available as they are, but have them convert from the same core Rust types.
- Because sorting is mutable API state for discrete histograms, the Python layer should expose both:
  - current ordered bins
  - stable bin identifiers on each discrete bin/snapshot entry
- For `evaluate_sample(s)` specifically:
  - returned histogram bundles must already have concrete labels
  - no process-aware unresolved label enums should escape through the Python API

### 8. Migrate observable runtime and tests

- Update `EventProcessingRuntime` and `Observables` to construct the new histogram core.
- Add Rust tests for:
  - discrete bin lookup
  - underflow/overflow behavior
  - labels in JSON
  - merge compatibility
  - rebin / contiguous-bin merge semantics
  - forced misbinning rejection
  - metadata-driven quantities with missing metadata
  - HwU compatibility against current discrete-emulation output
- Add Python API tests for:
  - standalone histogram creation/filling/merging
  - discrete snapshot inspection
  - JSON/HwU writing
  - event metadata exposure for `graph_name` / `graph_group_id` / `orientation_id` / `lmb_channel_id`
- Update architecture docs once the implementation shape is settled.

## Recommended Phasing

### Phase A: core design refactor

- Extract the histogram engine out of `HistogramObservable`.
- Introduce shared bin statistics and an explicit layout enum.
- Keep continuous behavior unchanged first.

### Phase B: discrete histogram support

- Add discrete layout/settings/snapshots/HwU writer.
- Wire discrete grouped filling and merge compatibility.
- Add separate `bin_id` vs rendered order.
- Implement no-op discrete `rebin(...)`.

### Phase C: discrete quantities and metadata plumbing

- Add `graph_group_id`, `orientation_id`, `lmb_channel_id`, and metadata-backed quantity definitions.
- Move count-like and cross-section-like quantities onto the discrete coordinate path.
- Add event metadata fields for the same identifiers.
- Add observable-local selections.
- Add automatic label and discrete-range resolution.

### Phase D: pyo3 exposure

- Wrap the standalone Rust histogram API.
- Add Python tests for direct histogram manipulation.

### Phase E: migration and cleanup

- Update current discrete-like builtins and tests.
- Update docs and examples.

## Remaining Design Points To Settle During Implementation

### 1. Exact settings shape for tagged histogram definitions

- Decided at high level: tagged continuous/discrete JSON.
- Still to settle in code:
  - whether the tag lives on `HistogramSettings` directly
  - or whether `ObservableSettings` flattens a tagged histogram-definition enum

### 2. Exact naming of the new quantity and metadata field

- User intent is clear: this is graph-group id, not graph id.
- Still to settle in code:
  - `graph_group_id`
  - `graphs_group_id`
  - another exact identifier
- Recommendation: use singular `graph_group_id` consistently.

### 3. Exact TOML/JSON shape for observable-local selections

- High-level direction is now clear.
- Still to settle in code:
  - exact tagged enum spelling
  - whether the inline quantity-based variant is literally `Selector(SelectorSettings)` or a sibling shape with renamed fields
- Recommendation:
  - keep it tagged and structurally close to existing selector definitions
  - keep it inline, not by name reference, so observables remain self-contained

### 4. Exact sort-strategy set for discrete histograms

- User requested:
  - ascending bin id by default
  - largest contribution first, both signed and absolute
- Still to settle:
  - whether sorting uses current average or raw sum
  - tie-breaking rules
- Recommended direction:
  - sort by current central value / signed contribution
  - or its absolute value
  - break ties by ascending `bin_id`

### 5. Scope of the public/internal metadata split

- Recommendation is to keep the split.
- Still to settle only if scope reduction is desired:
  - whether the split should apply to all serialized Rust evaluation outputs as well as Python
  - or whether only Python-facing conversion should expose `graph_name`

### 6. Cross-section default label policy

- User gave `total cross-section` as the intended kind of label.
- Still to settle:
  - whether this should be auto-filled by the quantity/histogram constructor
  - or only when requested through `bin_labels = Custom([...])`

## Immediate Next Step

- Pending review of the observable-local selections shape, the first implementation slice should be:
  - histogram-core extraction
  - tagged histogram settings refactor
  - shared event-predicate compilation usable by both global selectors and observable selections
  - internal sample-context plumbing for `graph_group_id` / `orientation_id` / `lmb_channel_id`
- While doing that, introduce the discrete sort model and the new stable `bin_id` representation immediately, because later retrofitting that separation would be expensive.

## 2026-04-04 Debugging Status

### LO cross-section regression after the discrete-histogram work

- User reported that
  - `./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate generate integrate -c "quit -o"`
    returned essentially zero on the real phase
  - while the card targets are `1.563e-3` and `5.329e-4`
- First investigation result:
  - this is not caused by observable filling, selector gating, or the new discrete histogram machinery
  - the exact same phase placement already existed on the parent revision right before the discrete-histogram commit
- Reproduced on the parent revision:
  - first iteration already gave `re ~ 1e-23`
  - and `im ~ -1.56e-3` / `-5.22e-4`
- Conclusion:
  - the physical LO result was already being accumulated into the imaginary component before the discrete-histogram branch

### Root cause found

- The LO path has threshold subtraction disabled, so the returned cross section is driven by the bare cut contribution.
- In `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`, the bare contribution assembled from `pass_two_result` had lost the subset-dependent complex phase.
- Restoring that phase at the bare contribution point fixes the LO phase placement:
  - from effectively pure imaginary
  - to effectively pure real

### Implemented fix

- Reintroduced the subset-dependent factor
  - `Complex::new_im(momentum_sample.one()).pow(num_esurfaces as u64)`
  - when converting `pass_two_result` into `bare_contribution`
- This is intentionally local to the cross-section bare-cut assembly, not to the general integrator or histogram accumulation logic.

### Verification

- Reduced check with `20k` samples and `--show-phase both`:
  - `epem_a_tth@LO`: `re ~ +1.59e-3`, `im ~ 0`
  - `epem_a_tth@LO_semi_inclusive`: `re ~ +5.17e-4`, `im ~ 0`
- Exact user command after the fix:
  - `epem_a_tth@LO`: `+1.547(21)e-3`
  - `epem_a_tth@LO_semi_inclusive`: `+5.330(72)e-4`
  - both within about `1σ` of the card targets

### Remaining follow-up

- The fix is verified on the LO `epem_a_ttxh` example and on the exact user command.
- Broader cross-section coverage outside this LO example should still be checked separately, but the immediate regression is resolved and is not related to the discrete histogram infrastructure itself.

## 2026-04-04 Observable Output Format Update

### User-requested behavior

- The user requested that
  - `default_runtime_settings.integrator.observables_output.format` accept a list of formats
  - e.g. `["hwu"]` or `["hwu", "json"]`
- The required semantics are:
  - the scalar form is no longer supported
  - when multiple formats are emitted, the final integration summary table must list all emitted observable files

### Implemented changes

- `ObservablesOutputSettings.format` now stores a list of observable file formats instead of a single one.
- Serde now accepts only the list TOML form:
  - `format = ["hwu"]`
  - `format = ["hwu", "json"]`
- The resolved output-format list is deduplicated and filters out `none`.
- The integration output writer now emits one latest observable file per configured format.
- The final "Integration results emitted" table now carries one row per emitted observable file and appends the format to the row label when multiple outputs are present for a slot.
- Checked-in TOML examples and fixtures that still used the scalar form were migrated to list syntax.

### Verification

- Rust settings tests:
  - single-entry lists parse and round-trip correctly
  - multi-format list parsing and round-tripping work
- Rust summary rendering test:
  - multiple emitted observable files for one slot render as separate rows
- End-to-end verification on
  - `./gammaloop --dev-optim --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate generate integrate -c "quit -o"`
  - with `format = ["hwu", "json"]` now reports both:
    - `integrands/epem_a_tth@LO/observables_final.hwu`
    - `integrands/epem_a_tth@LO/observables_final.json`
- Additional targeted integration tests in `tests/tests/test_runs/differential.rs` were updated to list syntax, but the full `test_runs` target currently does not compile because of an unrelated pre-existing missing field in `tests/tests/test_runs/events.rs`.

## 2026-04-04 Read-Only State and Run-History Follow-Up

### User-requested behavior

- When launching the CLI with `--read-only-state`, the CLI must never write into the active state folder.
- If a command would force such a write, it must fail cleanly and say that the operation cannot take place in `read-only-state` mode.
- In that mode, the default integration workspace must no longer live under the state folder:
  - use `./integration_workspace_<state_name>` when the state has a name
  - use `./integration_workspace` otherwise
- The persisted `run.toml` history must faithfully record all commands that actually executed, including:
  - boot-card `commands`
  - commands run through `run blockA blockB`
  - inline `-c "..."` commands
- The Python API must expose live in-memory state through:
  - `get_run_history()`
  - `get_global_settings()`
  - `get_active_command_blocks()`

### Implemented changes

- Added a shared read-only-state write guard in `crates/gammaloop-api/src/lib.rs`:
  - lexical path normalization
  - containment checks for targets under the active state
  - a single error path that explicitly mentions `--read-only-state`
- Applied that guard to:
  - `--clean-state`
  - default and explicit integration workspace creation
  - `save state`
  - `save dot`
  - `save standalone`
  - evaluator compilation triggered by `generate` when that compilation would write under the active state
- Changed the default read-only integration workspace location to:
  - `./integration_workspace_<state_name>`
- Flattened run-history recording in `crates/gammaloop-api/src/session.rs` so history captures the commands that actually ran instead of the wrapper `run ...` command.
- Added session-side TOML/string getters for the active in-memory run history and settings, plus command-block extraction.
- Exposed those through the Python API:
  - `GammaLoopAPI.get_run_history()`
  - `GammaLoopAPI.get_global_settings()`
  - `GammaLoopAPI.get_active_command_blocks()`

### Verification

- Rust-side checks:
  - `cargo fmt`
  - `cargo check`
  - `cargo test -p gammaloop-api clean_state_rejects_read_only_state_mode -- --nocapture`
  - `cargo test -p gammaloop-api one_shot_run_history_persists_all_executed_commands -- --nocapture`
  - `cargo test -p gammaloop-api read_only_state_rejects_workspace_inside_active_state -- --nocapture`
  - `cargo test -p gammaloop-api generate_rejects_compilation_into_active_state_in_read_only_mode -- --nocapture`
  - `cargo test -p gammaloop-api save_state_rejects_paths_inside_active_state_in_read_only_mode -- --nocapture`
  - `cargo test -p gammaloop-api save_dot_rejects_default_target_in_read_only_mode -- --nocapture`
  - `cargo test -p gammaloop-api save_standalone_rejects_paths_inside_active_state_in_read_only_mode -- --nocapture`
  - `cargo test -p gammaloop-api read_only_state_uses_cwd_workspace_default -- --nocapture`
  - `cargo test -p gammaloop-integration-tests --test test_cli cli_stateful_workflow_behaviors -- --nocapture`
- Python API validation:
  - refreshed the editable extension with `just build-api`
  - confirmed the live module exposes all three getters from `.venv`
  - `PYTHON=/Users/vjhirsch/Documents/Work/gammaloop_main/.venv/bin/python cargo test -p gammaloop-integration-tests --features python-api-tests --test test_python_api python_session_getters_return_live_toml_and_command_blocks -- --nocapture`

### Notes

- The Python rebuild issue turned out to be stale `maturin`/`cargo` processes from earlier interrupted attempts holding the target lock. After clearing those, the standard `just build-api` workflow completed normally and installed the updated extension into `.venv`.

## 2026-04-05 CLI Command Display Truncation

### User-requested behavior

- When the CLI logs `Running command <cmd>`, multiline commands longer than 5 lines should be truncated in the display.
- The displayed form should keep the first 5 lines and replace the remainder with a final `[...]` line.

### Implemented changes

- Added truncation logic in `crates/gammaloop-api/src/session.rs` inside the centralized command-display helper.
- The truncation applies only to multiline commands with more than 5 lines.
- Single-line commands and multiline commands up to 5 lines are displayed unchanged.

### Verification

- `cargo test -p gammaloop-api display_command_ -- --nocapture`
- `cargo fmt`
- `cargo check -p gammaloop-api`

## 2026-04-05 PR Review Follow-Up

### Review feedback addressed

- Validated discrete histogram ranges explicitly before allocation on both the Rust and Python histogram-construction paths.
- Fixed graph-group-context resolution for auto orientation/LMB domains and labels so it:
  - resolves by quantity type instead of literal selector quantity names
  - rejects negative singleton selector values
  - rejects out-of-range singleton graph/group ids when process info is available
  - avoids unchecked group indexing during label resolution
- Added a snake_case string helper for `DiscreteBinOrdering` and used it in Python snapshot serialization so ordering metadata round-trips correctly.
- Replaced the hot-path graph-group lookup-by-scan during graph evaluation with a cached graph-to-group mapping stored on the concrete process integrands.
- Cleaned the reviewed example-card nits:
  - removed the duplicated `remove processes` entry in `examples/cli/gg_hhh/1L/gg_hhh_1L.toml`
  - fixed the indentation issues that were introduced in the reviewed sections of `examples/cli/gg_hhh/1L/gg_hhh_1L.toml` and `examples/cli/aa_aa/1L/aa_aa.toml`
- Also changed the discrete selector matching branch to match by reference explicitly for consistency with the continuous selector path.

### Verification

- `cargo fmt`
- `cargo check`
- `cargo test -p gammaloop-api display_command_ -- --nocapture`

## 2026-04-05 Run Command Audit Assessment

### User-requested assessment

- The current saved `run.toml` records the unfolded commands executed by:
  - `./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate -c "quit -o"`
- The requested behavior is to keep the literal top-level command:
  - `run generate`
  rather than expanding the `generate` block contents into `commands = [...]`

### Assessment

- This is not actually specific to boot-card startup.
- The LO card does not define top-level `commands`; it only defines command blocks.
- The unfolding therefore comes from the generic `run ...` command-history behavior in `CliSession::execute_run(...)`, not from `apply_boot_run_history(...)`.

### Difficulty

- If the goal is only cosmetic or UX-oriented, and we accept a change in the audit model, then the implementation is moderate and localized:
  - introduce a distinction between:
    - externally entered top-level `run ...` commands
    - nested/internal `run ...` expansion
  - record the wrapper form for the top-level case
  - keep nested runs flattened
  - update the existing run-history tests accordingly
- If the goal is to keep the current “faithful replay/audit” property in the presence of later command-block redefinitions, then it becomes more invasive.
  - Today the flattening guarantees that replay does not depend on the later meaning of a block name.
  - Recording only `run generate` would make replay depend on the stored `command_blocks` definitions for `generate`.
  - That is acceptable only if we are comfortable with that semantic shift, or if we add a stronger snapshotting/versioning model for command-block definitions at execution time.

### Recommendation

- I would classify the simple behavior change as moderate, not hard.
- The tricky point is not implementation complexity but semantics:
  - whether we still want `run.toml` to be an execution transcript of concrete commands
  - or whether we want it to preserve higher-level `run <block>` invocations even though that makes replay depend more directly on block definitions
- `cargo test -p gammalooprs discrete_histogram_ -- --nocapture`
- `cargo test -p gammalooprs resolve_graph_group_context_uses_selector_quantity_type -- --nocapture`
- `cargo test -p gammalooprs resolve_graph_group_context_rejects_negative_singleton_values -- --nocapture`

## 2026-04-05 Run Command Audit Implementation

### User-approved direction

- The higher-level invocation form is acceptable.
- Command-block definitions are treated as stable enough that `run.toml` can preserve:
  - `run generate`
  instead of unfolding the block contents.
- The user also noted this was already the intended behavior for CLI-entered `run blockName...` commands, so the implementation should restore that wrapper-level recording model consistently.

### Implemented changes

- Changed `CliSession::execute_run(...)` in `crates/gammaloop-api/src/session.rs` so that:
  - child commands inside a `run ...` execution are always evaluated with `HistoryMode::Suppress`
  - the wrapper `run ...` command itself is what gets recorded when the outer history mode is `Record`
  - this also applies when the run plan exits through `quit -o` / `ControlFlow::Break(...)`
- Added persisted-history sanitization for wrapper `run ...` commands:
  - non-persistable inline commands are removed before the wrapper is saved
  - in practice this means `quit`, `start_commands_block`, and `finish_commands_block` never appear inside a persisted `run ... -c ...`
  - if sanitization removes the entire inline `-c` payload, the history entry falls back to the plain wrapper such as `run generate`
- This restores a higher-level audit script for:
  - direct CLI `run ...`
  - one-shot `... run ...`
  - boot-run-history replay of saved `run ...` entries
  - interactively defined blocks later executed through `run <name>`
- Nested runs are no longer flattened into their child commands in the saved history; the top-level recorded wrapper remains the single audit entry.

### Tests updated

- Updated the run-history expectations in:
  - `tests/tests/test_cli.rs`
  - `crates/gammaloop-api/src/lib.rs`
- Specifically, the following now expect wrapper-form history entries:
  - `run set_display set_runtime -c "..."`
  - `run outer`
  - `run boot_block`
  - `run demo`
  - `run demo -c 'quit -o'` now persists as `run demo`
  - one-shot `run cmdBlockA cmdBlockB -c "...; quit -o"` now persists without the trailing `quit`

### Verification

- `cargo fmt`
- `cargo test -p gammaloop-api one_shot_run_history_persists_all_executed_commands -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_cli cli_stateful_workflow_behaviors -- --nocapture`
- `cargo check`

## 2026-04-06 Amplitude Observable Assessment

### User request

- Assess how much work it would be to support observables in the amplitude evaluation stack as well.
- In particular:
  - discrete observables such as total contribution, per-graph, per-orientation, and per-LMB-channel should also make sense for amplitudes
  - `evaluate_sample(s)` should then be able to return `HistogramBundle` snapshots for amplitudes too
  - particle-derived quantities may be interpreted from the fixed final state of the amplitude process, even if they are effectively constant across samples

### Deep-dive conclusion

- This is smaller than a new histogram feature. The result types, observable bundles, Python wrappers, and process-level snapshot plumbing are already generic over both process-integrand kinds.
- The main missing piece is amplitude-side event production policy and selector handling, not bundle transport.

### What is already in place

- `ProcessIntegrand::observable_snapshot_bundle(...)`, `build_observable_snapshots_for_result(...)`, and the evaluation result types already support amplitudes and cross sections uniformly.
- `AmplitudeIntegrand::warm_up(...)` already constructs an `EventProcessingRuntime` with:
  - selectors
  - observables
  - histogram process info for graph/group/orientation/LMB labels
- Amplitude graph terms already know how to generate an event carrying:
  - incoming/outgoing PDGs
  - outgoing external momenta
  - `orientation_id`
  - `lmb_channel_id`
  - `lmb_channel_edge_ids`
- Generic process-level code already injects:
  - `graph_id`
  - `graph_group_id`
  into every returned event group for both amplitudes and cross sections.
- Observable quantity extraction is already generic:
  - metadata quantities (`graph_id`, `graph_group_id`, `orientation_id`, `lmb_channel_id`)
  - particle quantities
  - jet quantities
  - `cross_section` / total-weight style observable

### Main blockers identified

- `AmplitudeGraphTerm::evaluate(...)` currently ignores the `event_processing_runtime` argument entirely.
- It only generates an event when `settings.should_return_generated_events()` is true, i.e. when `general.generate_events` is explicitly enabled.
- That is inconsistent with the runtime policy encoded in `RuntimeSettings`:
  - `should_generate_events()` is true when active selectors or observables exist
  - `should_buffer_generated_events()` is true when events should be retained for outputs or observables
- Because of that mismatch:
  - amplitude observables do not naturally work in integration runs where observables exist but `general.generate_events = false`
  - active selectors are also not really enforced in the amplitude stack

### Practical implementation scope

- I would classify this as moderate and fairly localized.
- The likely implementation is:
  - mirror the cross-section event-generation policy inside `AmplitudeGraphTerm::evaluate(...)`
  - build amplitude events whenever selectors or observables require them
  - run selectors through `EventProcessingRuntime`
  - zero the graph contribution when an active selector rejects the event
  - only buffer identity-rotation events for returned event groups / observable accumulation
- No result-type redesign appears necessary.
- `evaluate_sample(s)` in the Python/Rust API already returns observable bundles generically; once amplitude event groups are produced under the right conditions, those APIs should start returning non-empty bundles automatically.

### Likely short-comings / caveats

- Particle and jet quantities for amplitudes are only as meaningful as the fixed external-state interpretation:
  - they will often be constant across samples
  - this matches the user’s proposed acceptance criterion, but is still semantically weaker than cross-section event observables
- `cross_section` as a quantity name is slightly misleading on amplitudes:
  - technically it would histogram the sampled amplitude contribution / integrand weight
  - functionally this is still useful, but the naming may deserve a later alias such as `integrand_weight`
- Active selector semantics on amplitudes need to be agreed implicitly as:
  - “accept/reject the entire graph contribution represented by that event”
  - not “accept/reject individual cuts”, because amplitudes have no cut-level event decomposition
- The precise evaluation path (`evaluate_sample_precise(s)`) currently builds per-sample observable snapshots from returned event groups, but unlike the standard Python `evaluate_sample(s)` path it does not force `general.generate_events = true`.
  - if precise amplitudes should also expose observables reliably, that path should be aligned too
- Rotated stability evaluations should continue to avoid double-counting observables:
  - the safest approach is to mirror the cross-section rule and only buffer identity-rotation events, while using rotated events only if selector checks truly need them

### Recommended implementation order

- 1. Make amplitude event-generation policy follow `should_generate_events()` / `should_buffer_generated_events()` instead of `should_return_generated_events()` alone.
- 2. Apply selectors inside `AmplitudeGraphTerm::evaluate(...)` using `EventProcessingRuntime`, mirroring the cross-section flow as closely as possible.
- 3. Add amplitude-focused tests:
  - active selectors reject amplitude contributions
  - metadata quantities (`graph_id`, `graph_group_id`, `orientation_id`, `lmb_channel_id`) fill discrete histograms
  - particle quantities work on the fixed external state
  - API `evaluate_sample(s)` returns non-empty observable bundles for amplitudes
- 4. Optionally align the precise evaluation API with the same forced-event behavior if that is desired.

### Overall assessment

- The user’s intuition is correct: this should not be too difficult.
- The observable infrastructure is already generic enough.
- The real work is mostly:
  - amplitude event buffering policy
  - amplitude selector application
  - tests and a small amount of semantic clarification

## 2026-04-06 Amplitude Observable Revised Plan

### User clarifications

- Amplitudes must incur zero event-generation overhead when:
  - `general.generate_events = false`
  - no selectors are active
  - no observables are configured
- If selectors or observables are present, `generate_events` becomes irrelevant for internal evaluation:
  - events must be generated internally
  - exactly as in the cross-section flow
- Code duplication between amplitude and cross-section event handling should be reduced rather than increased.
- The observable quantity currently named `CrossSection` should be renamed to `Integral`, with all TOML cards updated accordingly.
- The current `evaluate_sample(s)` behavior should be reassessed carefully:
  - users should still be able not to receive returned event groups
  - internal event generation for selectors/observables must remain possible without surfacing those events to the caller

### Revised deep-dive conclusion

- The runtime already contains the right conceptual split:
  - `should_generate_events()` for internal event creation
  - `should_buffer_generated_events()` for keeping accepted identity-rotation events around long enough for observables / caller-facing output
  - `should_return_generated_events()` for whether event groups survive all the way into returned API results
- Cross sections already implement this split correctly.
- Amplitudes currently do not:
  - they only generate events when `should_return_generated_events()` is true
  - they ignore `event_processing_runtime` inside graph evaluation
- The Python / inspect `evaluate_sample(s)` helper layer currently muddies the distinction by temporarily forcing `general.generate_events = true`, which forces both internal generation and returned event groups.

### Current control-flow status

- In the process-evaluation core:
  - results are passed through `process_evaluation_result_runtime(...)` before any possible result-side event discard
  - then `maybe_discard_generated_events_in_result(...)` drops returned event groups when `should_return_generated_events()` is false
- This means the low-level process runtime already supports:
  - internal event generation for selectors/observables
  - without necessarily surfacing event groups to the caller
- Cross-section graph evaluation uses that model already.
- Amplitude graph evaluation is the outlier.
- Separately, `EvaluateSamples` in the API layer currently has `force_generate_events: bool`, and when true it mutates `general.generate_events = true` on the integrand before evaluation.
  - that is why Python `evaluate_sample(s)` and CLI `inspect` currently always get returned events
  - this is too blunt if we want internal generation without mandatory surfacing

### Recommended implementation approach

- 1. Align amplitude graph evaluation with cross-section event policy.
  - Generate events when `settings.should_generate_events()` is true.
  - Buffer only when `settings.should_buffer_generated_events()` is true.
  - Run selectors via `EventProcessingRuntime`.
  - Use rotated events only for selector checks when needed, not for observable/event buffering.
- 2. Refactor the shared event-policy logic into common helpers in `integrands/process/mod.rs`.
  - Do not duplicate the cross-section event-decision flow into amplitudes.
  - Share the logic for:
    - deciding whether an event must be built
    - deciding whether it is buffered
    - running selectors
    - updating generated/accepted counters
    - dropping non-buffered / selector-rejected events
  - Keep only the graph-specific event-construction details in each integrand kind.
- 3. Separate “internal generation” from “returned events” in the API helper layer.
  - Replace or refine `force_generate_events` so it no longer means “mutate `general.generate_events = true`”.
  - The API request should distinguish:
    - force internal event generation for inspection / selectors / observables
    - whether returned event groups should be surfaced to the caller
  - This lets:
    - amplitudes and cross sections behave identically
    - Python `evaluate_sample(s)` optionally avoid returning event groups
    - internal observable/selector logic still function
- 4. Rename the observable quantity `CrossSection` to `Integral`.
  - Update:
    - `QuantitySettings`
    - completion/schema/catalog code
    - serialization names
    - tests
    - example cards
  - No backward-compatibility alias is needed.
- 5. Add amplitude observable coverage.
  - Metadata discrete histograms:
    - `graph_id`
    - `graph_group_id`
    - `orientation_id`
    - `lmb_channel_id`
  - Constant external-state particle quantities.
  - Active-selector rejection behavior.
  - API `evaluate_sample(s)` observable bundles for amplitudes.

### Code structure recommendation

- The most valuable deduplication target is not the numerical integration logic itself.
- It is the shared event-handling policy around a freshly generated candidate event.
- A good structure would likely factor out a helper that accepts:
  - whether the current rotation is identity
  - runtime settings
  - optional event-processing runtime
  - a closure or already-built event
- and returns:
  - selector acceptance
  - maybe-buffered event
  - event-processing timing
  - generated/accepted count deltas
- Cross-section can call it per cut-event.
- Amplitude can call it per graph-contribution event.

### Expected shortcomings / caveats

- Particle/jet observables on amplitudes may be sample-constant today because amplitude external kinematics are typically fixed, but the implementation should not encode that assumption.
- In particular, we should not introduce any cross-sample caching or special-case shortcut for amplitude particle/jet observables based on current constant externals.
- The observable path should continue to read event kinematics directly so it remains correct if amplitudes later support sample-dependent external momenta.
- `Integral` is more correct than `CrossSection`, but some display text and docs will need to be adjusted so users understand it now covers both cross sections and amplitudes.
- If exact semantic parity is desired, the precise evaluation API should be reviewed in the same patch:
  - it currently builds observable snapshots from event groups directly
  - and does not use the same top-level `force_generate_events` hook as the standard API path

### Recommended rollout order

- 1. Rename `CrossSection` quantity to `Integral`.
- 2. Refactor shared event policy into common helpers.
- 3. Make amplitude evaluation use the shared policy.
- 4. Refine `EvaluateSamples` / Python inspection flow to separate internal generation from returned event surfacing.
- 5. Add amplitude observable and selector tests.

### Overall assessment

- Still moderate.
- The main extra work introduced by the user clarifications is not amplitude observables themselves.
- It is the API-layer cleanup needed to stop conflating:
  - internal event generation
  - buffering for observables/selectors
  - surfacing event groups to the caller

### Final revision after external-kinematics clarification

- Keep particle, jet, and continuous observables fully supported for amplitudes even if they are mostly trivial with today’s constant externals.
- Do not add any amplitude-specific caching of derived observable inputs across samples.
- The implementation should stay event-driven:
  - build the event when runtime policy says it is needed
  - let the observable/selector machinery read directly from that event
  - avoid baking in the assumption that amplitude externals are constant across samples
- This keeps the design future-proof for the existing `constant` vs non-constant external-momenta model without adding extra complexity now.

## 2026-04-06 Amplitude Observable Implementation

### What was implemented

- Added a shared prepared-event helper in `crates/gammalooprs/src/integrands/process/mod.rs`:
  - `prepare_buffered_event(...)`
  - `PreparedBufferedEvent<T>`
- This helper centralizes the common event policy for process integrands:
  - whether an event must be built
  - selector processing
  - whether an accepted event is buffered
  - generated / accepted event counters
  - event-processing timing
- Cross-section evaluation was refactored to use that helper instead of keeping a local duplicate event-processing path.
- Amplitude evaluation was refactored to use the same helper, so amplitudes now honor selectors and observables with the same internal event-generation policy as cross sections.

### Resulting event policy

- When `general.generate_events = false` and there are no active selectors and no observables:
  - amplitudes stay on the event-less fast path
  - there is no additional event-generation overhead
- When selectors or observables are present:
  - amplitudes now generate internal events even if `general.generate_events = false`
  - selectors are evaluated through the shared event-processing runtime
  - rejected amplitude samples have their whole graph contribution zeroed
  - accepted identity-rotation events can still feed observables even when events are not surfaced to the caller
- Returned event groups remain controlled separately from internal generation.

### API control-flow cleanup

- Replaced the old `force_generate_events` knob in the API helper layer with:
  - `return_generated_events: Option<bool>`
- This was applied to:
  - `EvaluateSamples`
  - `EvaluateSamplesPrecise`
  - the Python `evaluate_sample(...)` and `evaluate_samples(...)` bindings via a `return_events` argument
- `inspect` now explicitly requests returned events, instead of relying on the old blunt forcing behavior.
- Internal event generation is therefore now governed by runtime need, while caller-visible returned events are separately controllable.

### Observable quantity rename

- Renamed the observable quantity `CrossSection` to `Integral` throughout:
  - observable settings / definitions
  - completions and command templates
  - tests
  - example cards
- Updated the example cards in:
  - `examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml`
  - `examples/cli/epem_a_ttxh/NLO/epem_a_tth_NLO.toml`
  - `examples/cli/epem_a_ttxh/LO/OLD_no_grouping_individual_top_antitop_epem_a_tth_LO.toml`

### Test and fixture updates

- Migrated test-side `jet_count` observables to the new discrete histogram form.
- Updated the amplitude CLI fixture helper to the current sampling syntax:
  - `graphs = "monte_carlo"`
  - `orientations = "monte_carlo"`
  - `lmb_multichanneling = true`
  - `lmb_channels = "monte_carlo"`
- Added amplitude regression coverage for:
  - selector-driven internal event generation without surfacing events
  - observable accumulation without returned events
  - per-call override of whether returned events are surfaced

### Documentation updates

- Updated `docs/architecture/architecture-current.md` to document the shared event policy and the split between:
  - internal generation
  - buffering
  - returned event surfacing
- Updated `IGNORE/differential_lu.md` to reflect that amplitudes now use the same prepared-event path as cross sections for selector and observable handling.

### Validation

- `cargo fmt`
- `cargo check`
- `cargo check -p gammaloop-api --features ufo_support,python_api`
- `cargo test -p gammalooprs observable_settings_infer_discrete_histogram_from_domain -- --nocapture`
- `cargo test -p gammaloop-api completion_offers_quantity_kinds_for_process_add -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_runs amplitude_events_surface_threshold_counterterms_and_reproduce_weight -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_runs amplitude_selectors_generate_internal_events_without_surfacing_them -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_runs amplitude_observables_accumulate_without_returning_events -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_evaluation_api lu_rust_evaluate_samples_can_override_returned_events_without_forcing_internal_generation -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_evaluation_api lu_rust_evaluate_samples_respect_event_generation_and_observables -- --nocapture`
- `cargo test -p gammaloop-integration-tests --test test_runs inspect_x_space_reports_missing_discrete_dimensions_cleanly -- --nocapture`
- `cargo clippy --all-targets`

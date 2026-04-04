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

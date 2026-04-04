# Discrete Histograms Human and Codex Interactions

## 2026-04-04 User

The current project gammaloop has support for observables.
Some of them, like the cross-section and jet count is discrete in nature, and for now emulated through the normal continuous histograms with x-axis boundaries given by continuous values.

I would like to change this and make discrete-value histograms first-citizen, but still using a similar API as existing histograms.
Moreover, I would like the whole histograms infrastructure (creating them, filling them in with correlated grouped weights, merging bings, merging histograms, etc...) to be exposed using pyo3 through our python API.
This discrete variant of histograms must still write out to the HwU format in exactly the same way, but of course the json output may be slightly different.

What I want to have better handled for this new discrete variant of histograms is the following:

- There is no nothing of bin width, each bin by definiting covers one single unit. It may however still have a start and end, with spill-over recorders (below and above).
- I want to have the possibility to add names to the bin (i.e. "String identifiers" that can be used later for rendering.
- Misbinning treatments must always be disabled for such histograms.
- Regarding filling, ideally the x-value should be supplied through an isize. However, for structural/API reason, if this needs to be still done through floats that is acceptable too.
- I will want to add hardcoded versions of "quantities" that are discrete in nature (similar to current jet count), such as "orientation_id", "lmb_channel_id" and "graph_id", that can be used as input for such discrete histograms.
Events metadata should also be completemented with such discrete usize quantities, but optional because dependening on the sampling strategies those do not make sense (and therefore the quantity should be none, and histograms not filled). I want you to investigate the details of how such information (orientation, graph and lmb) ids can be forwarded efficiently and elegantly down to these quantities and events metada, and how can filling of histograms can be done elegantly in such cases where these discrete quantities may or may not be provided (empty list if not available I guess)?

Anyway formulate a first plan for all the points above, with identification of the tricky points that require common review together. We'll review the resulting plan together.
During the whole planning and implementation procedure, keep track a details status and objectives in `discrete_histograms.md` and also collect all our interactions (my prompts and your answers, but not your intermediate investigations) in `discrete_histograms_human_and_codex_interactions.md`
And keep those files up to date at all times.

## 2026-04-04 Codex

First planning response:

- Refactor the histogram implementation so the accumulation/fill/merge/snapshot/HwU logic is a standalone histogram engine instead of being embedded inside `HistogramObservable`. This is the prerequisite for both first-class discrete histograms and a clean pyo3 API.
- Introduce an explicit histogram layout split between continuous and discrete axes. The discrete layout should carry `start`, `end`, optional per-bin labels, underflow/overflow bins, and no misbinning support.
- Introduce explicit observable coordinates so discrete quantities use a real integer path internally:
  - continuous entries stay continuous
  - discrete entries use `isize`
- Add discrete-native quantities for `graph_id`, `orientation_id`, and `lmb_channel_id`, and move naturally discrete existing quantities such as counts and cross-section onto the discrete coordinate path.
- Extend event metadata with optional `orientation_id` and `lmb_channel_id`, while keeping existing detailed `lmb_channel_edge_ids`. When optional metadata is absent, the quantity should emit no entry so the histogram is simply not filled.
- Keep HwU output numerically identical to the current continuous-emulation convention by writing integer bin `k` as `[k, k + 1)` and a single cross-section bin as `[0, 1)`.
- Expose the standalone histogram core through pyo3 with constructors and methods for grouped filling, merge, update/finalize, snapshot/restore, rebin or contiguous-bin merge, and JSON/HwU writing.

Points requiring common review before implementation:

- Whether to make a clean breaking change for old count/cross-section run cards, or keep a temporary compatibility layer.
- The JSON shape for discrete histograms.
- The semantics of rebinning or bin merging when discrete bins carry labels.
- Confirmation that `[k, k + 1)` is the intended exact HwU convention for discrete bins.
- Whether `lmb_channel_id` should be the compact multi-channel index while keeping edge-id lists separately.
- Confirmation that missing optional metadata should mean "emit no entry" rather than using sentinel bins.
- Whether cross-section should get a default visible label such as `total`.

Detailed status and objectives are tracked in `discrete_histograms.md`.

## 2026-04-04 User

1. Clean break, existing jet count and cross-section should be migrated to this new model, as well as all toml cards defining them that you see.

2. Make sure to adjust CLI completions for this new type of quantities and histograms when/if need be.

3. In JSON I think I prefer tagged discrete / continuous shape.

4. rebinning discrete histograms leave them unchanged.

5. You can lock these conventions for HwU.

6.`lmb_channel_id` means the same as the value of the discrete bin supplied when doing lmb multi-channeling sampling. "lmb_channel_edge_ids" is a completely separate thing (which will be used to label discrete bins, but see later point for this).

7. No metadata available (because of sampling strategy typically) means no fill. Whether you go through empty list or a None value at intermediate steps, is up to you (I lean towards empty list).
For the metadata of the event however, I would like None to be there if not available.

8. Regarding bin-labels in general, I would like that the definition of discrete histograms have the possibility to supply descriptive labels for each bins (through a list that will be ordered in the same way as the discrete bins within the interval [min:max]. So for cross-section I could simply think of writing there something like "total cross-section".
However, for orientations, lmb_channels and graphs, I obviously would like to avoid having to supply them manually. So I was thinking of the field "bin_labels" of such distrete histograms to take the variant `Custom(["bin_label_1", "bin_label_2", ...])`,`Orientation`,`GraphsGroupMasterName` or `lmb_channel_edge_ids` in which case these bin_labels strategies would automatically create the label for the corresponding orientations (using the uncolored version of the conventions +-0 notation for each edge that you can find elsewhere in the code e.g. in the display command), the name of the group master graph for this group ID, or a string indicating the sorted list of lmb_edge_ids.
Notice again that before when I mention `graph_id`, I really meant `graphs_group_id`, i.e. as always mapping exactly the meaning of the discrete sampling dimension value when enabled (otherwise not available, as discussed earlier).

9. Another important aspect of the implementation of the discrete histograms is that I want to be able to sort the bins according to diffferent strategies (i.e. largest contribution first (in abs value, or directly), or "original" i.e. following ascending bin id, etc.... Such strategies should be possible to specify during definition of the histograms and also be changed at any point later as part of the API.
This means that the bin order in the list does *not* necessarily match the bin ID, so that it needs to be kept separate because obviously during filling and merging etc... one must keep the information of the bin id independently of the position of the bin in the list (different than the bin labels, which could in principle be identical even, those are just meant fro rendering).

## 2026-04-04 Codex

Refined planning response:

- I will treat this as a clean migration, including the existing `jet_count` / `cross_section` paths and the TOML cards in the repository that define them.
- I will update the settings/completion surface as part of the type changes, because completion is driven from the serialized default settings plus schema.
- The histogram model should now be explicitly tagged between continuous and discrete definitions in JSON and settings.
- Discrete `rebin(...)` will be implemented as a no-op, returning the histogram unchanged.
- The new metadata-backed discrete quantity is `graph_group_id`, not `graph_id`, to match the discrete sampling dimension semantics. Event metadata should therefore keep concrete `graph_id` and gain optional `graph_group_id`, `orientation_id`, and `lmb_channel_id`.
- Missing metadata will remain `None` on events, while quantity extraction will simply emit no entries and therefore not fill the histogram.
- Discrete histograms need label strategies:
  - `Custom([...])`
  - `Orientation`
  - `GraphsGroupMasterName`
  - `LmbChannelEdgeIds`
- For automatic labels, I will reuse current naming conventions already present in the code:
  - uncolored `+-0` strings for orientations
  - master graph names for graph groups
  - sorted edge-id list strings for LMB channels
- Discrete histograms also need a strict separation between:
  - stable `bin_id` used for filling and merging
  - current rendered order used in JSON and Python APIs
- I will therefore design discrete histograms around stable bin ids plus a configurable sort strategy such as:
  - original ascending id
  - largest contribution first
  - largest absolute contribution first

Remaining implementation-time design choices are narrower now:

- the exact Rust/settings type shape for the tagged histogram enum
- the exact spelling of `graph_group_id`
- how automatic label strategies are supported for standalone pyo3 histograms that may not have process metadata attached
- the exact contribution metric and tie-breaks used by sort strategies
- whether `cross_section` gets an automatic default label or only an explicitly configured one

Detailed status and the refined plan are tracked in `discrete_histograms.md`.

## 2026-04-04 Codex

Additional clarification after refining the plan:

- I am locking the discrete histogram range assumption to an inclusive integer interval `[min:max]`.
- That means the in-range discrete bin ids are exactly `min, min + 1, ..., max`.
- This matches the requested label-list behavior and stays compatible with HwU output by writing each integer bin `k` as `[k, k + 1)`.

## 2026-04-04 User

A couple clarification.

- For the event metadata I would like the graph_id to actuall be graph_name directly, but the graph_group_id can rename a usize graph_group_id indeed.
- For the Python API it does not make sense to allow graph labels option to be the auto-label variants, and they must have been concretised at the point they are used to be returned as part of the histogramsbundle handed back to the caller of evaluate_sample(s).
- Make sure that when calling `evaluate_sample(s)` if the information regarding orientation, graph_group, or lmb_channel_id is being supplied in any way (I guess most likely through discrete dimension values, but it may be something else too maybe, check!) then it must find itself in the event metadata and also in the corresponding discrete histograms to be filled.
- Add an auto-label varian "BinID" which assigns the label `#<bin_id>` to each of such discrete bin.
- Keep in mind that whenever you are within the evaluation stack with event generation or histogram filling, this is a very efficiency sensitive part of gammaloop code and you must make sure to implement things as efficiently as possible and with as little allocation as possible.

Revise your plan based on the above, do a deep re-assessment and come back with an updated plan and possibly additional questions.

## 2026-04-04 Codex

Updated planning response after deeper reassessment:

- I will treat internal event metadata and returned/public event metadata as separate layers:
  - internal event metadata should stay compact and numeric so the hot evaluation and histogram-filling path can avoid repeated string allocation and lookup
  - returned event metadata should expose `graph_name` plus optional `graph_group_id`, `orientation_id`, `lmb_channel_id`, and `lmb_channel_edge_ids`
- I will normalize graph-group, orientation, and channel selections after they have already converged into `GammaLoopSample` / `DiscreteGraphSample`, rather than trying to special-case each frontend input route. That guarantees x-space `discrete_dims`, momentum-space `discrete_dims`, explicit momentum `graph_name`, and explicit momentum `orientation` all populate metadata and fill the corresponding discrete histograms consistently.
- The histogram core refactor remains the first implementation step, but it now needs an explicit split between:
  - tagged continuous vs discrete histogram definitions
  - stable discrete `bin_id` storage
  - mutable rendered order through discrete sort strategies
- Discrete label strategies should now be modeled as a settings-time concern only:
  - `Custom([...])`
  - `BinID`
  - `GraphsGroupMasterName`
  - `Orientation`
  - `LmbChannelEdgeIds`
  but snapshots returned from `evaluate_sample(s)` and the Python-facing histogram bundle must only contain concrete resolved labels.
- Process-aware labels should therefore be resolved once, preferably during observable compilation/runtime initialization when the integrand is known, not during per-event filling.
- Standalone pyo3 histogram construction should expose only concrete labels; even `BinID` can be concretized immediately there.
- The main new tricky point is that `orientation_id` and `lmb_channel_id` are local to a graph group, not globally meaningful across all graph groups. Their automatic label strategies therefore become ambiguous if one histogram can receive contributions from multiple graph groups.

Main open review question from this reassessment:

- What should happen when an `orientation_id` or `lmb_channel_id` histogram can receive contributions from multiple graph groups, where the same bin id means different things?
  - My current recommendation is to allow the quantities themselves, but only allow `Orientation` / `LmbChannelEdgeIds` automatic labels when the histogram is effectively tied to one graph group, or when GammaLoop can prove the mapping is identical across all contributing graph groups.

Detailed status and the updated reassessment are tracked in `discrete_histograms.md`.

## 2026-04-04 User

5. How complicated is the split of public and internal event metadata? If it is too complicated then it may not be worth it.

7. Default should be "ascending_bin_id" (so static). But allow for an API "change_bin_ordering" which changes this variant and consequently re-order the bin (though filling / merging etc.. should always work)>

Regarding your point about the graph-group-local nature of the Orientation and LmbChannelEdgeIds this is indeed very true and we must address this carefully.
I would like to solve this by allowing *all* observables to now revieve another field called "selections" (defaulting to empty). which contains a list of "Selection" variants allowing the user to add requirements for the observable to be filled in by such event.
The admissible Selection variant can be the definition of a selector similar to the existing Selectors that can be added to the run details. Or it can be "GraphGroupIDEqual(<graph_group_id_selected>" or "OrientationIDEqual(<orientation_id_selected>)" or "LMBChannelIDEqual(<lmb_channel_id_selected>).
In turns, such selectors will allow to the relevant auto-labels variants to resolve and assign proper variants using the information from the `GraphGroupIDEqual` (if absent then crash elegantly).
Same for the bin integer range, instead of allowing the user to just add [min, max] allow it to supply a variant that resolves to maximum number of orientations, lmb_channel_ids or graph groups, i.e. `[0, NMaxOrientations]`. These variants would again not resolve if the `GraphGroupIDEqual` is absent and there are more than one.

I know this last point above is quite tricky, so do another deep pass of the plan and what you concluded, and come back with an update for a review together.

## 2026-04-04 Codex

Updated planning response after the observable-selection design pass:

- The public/internal event-metadata split is moderate in complexity, not tiny but still contained:
  - `EvaluationResultOutput` currently serializes `EventGroupList` directly
  - `EventGroupList` uses the same `CutInfo` type as the internal hot path
  - so a true split requires dedicated output event/cut-info structs plus conversion in `EvaluationResult::into_output(...)` and the precise equivalents
  - it does not require changing histogram filling logic itself
- My current recommendation is still to keep that split, because otherwise `graph_name` would need to be allocated and carried per event in the hot event-generation path.
- I will switch the default discrete ordering to `ascending_bin_id`, and treat `change_bin_ordering(...)` as an API operation that only changes rendered order while fill/merge compatibility remains keyed by stable `bin_id`.
- Your new `selections` proposal solves the graph-group-local label problem cleanly. I now recommend:
  - add `selections: []` to all observables
  - support an inline selector-definition variant reusing the existing selector logic
  - support direct metadata predicates `GraphGroupIDEqual`, `OrientationIDEqual`, and `LMBChannelIDEqual`
  - AND them together before the observable extracts entries from an event
- I also recommend using those observable-local selections as the resolution context for:
  - `Orientation` and `LmbChannelEdgeIds` auto labels
  - auto discrete ranges such as graph-group ids, orientation ids, and channel ids
- For the auto ranges, because discrete ranges are inclusive, the concrete resolution should be:
  - graph groups: `[0, graph_group_count - 1]`
  - orientations for a selected group: `[0, n_orientations(group) - 1]`
  - channels for a selected group: `[0, n_channels(group) - 1]`
- Resolution rule for group-local auto labels/ranges:
  - use `GraphGroupIDEqual(...)` if present
  - otherwise, if the integrand has exactly one graph group, use that implicitly
  - otherwise fail early with a clear error

Detailed status and the updated plan are tracked in `discrete_histograms.md`.

## 2026-04-04 User

1. Actually no, not inline definitoins, I think simplest would be to allow selectors to have a flag deciding whether they are "active" or "inactive" as far as integrand runtime event selection goes (i.e. their current role).
This way, it allows one to register them but use them by referring to their name only in the context of the observables selection, not the main event generation.

2. I prefer to no have the split of internal and output metada for the Event. So instead simplify the metadata to be kept numerical (i.e. no graph_group_master_name or graph_name, but the correspondings "ids" only).
This avoids the complexity of separating internal and output-facing metadata. But please comment on whether this is enough to be able to keep internal and output-facing data the same without loss of performance.

3. Regarding the range, you are right, the -1 in [0, n_orientations - 1] is necessary, but it's a bit awkward. My objective is just to allow for the user not having to manually fetch these ranges for the orientations, graphs and lmb_edge_ids discrete histograms. If you can come up with anything more elegant to solve that problem, then please suggest it.

Review the plan thoroughly be re-reading all relevant code and assessing the validity of what has currently been decided as well as assessing the points above.
Come back then with the updated plan and we'll review it.

## 2026-04-04 Codex

Updated planning response after the selector-reference and shared-metadata pass:

- Switch observable-local selections to named selector references rather than inline selector definitions.
- Add `active: bool` to selectors so they may stay registered and reusable while being ignored by the runtime event-filtering path.
- Drop the internal/output event-metadata split and keep one shared numeric `CutInfo`, which is sufficient for performance because it avoids hot-path string allocation and only adds small fixed-size integer fields.
- Hide the awkward explicit `[0, n - 1]` user spelling behind discrete-domain variants such as `GraphGroupIds`, `OrientationIds`, `LmbChannelIds`, and `SingleBin`.
- Keep group-local auto-label and auto-domain resolution tied to singleton graph-group context, resolved from referenced selectors or from the single-group-integrand case.

Detailed status and the revised plan are tracked in `discrete_histograms.md`.

## 2026-04-04 User

I would like you to keep graph_id and graph_group_master_id as two separate and available discrete quantities for the discrete histograms. The corresponding autolabel and autorange variants for these two separate quantities must exist too then. Also make sure that selectors can also support discrete quantities in input and allow these discrete quantities to be within some user specified inclusive range.

Ask me questions about what is still unclear given the updated plan you'll communicate.

## 2026-04-04 Codex

Updated planning response after reviewing the additional quantity and selector requirements:

- Keep `graph_id` as a discrete quantity in addition to `graph_group_id`.
- Add selector support for discrete quantities through an inclusive `DiscreteRange { min, max }` variant.
- Keep observable `selections` as named selector references only.
- Add auto-domain variants for `GraphIds` and `GraphGroupIds`.
- Add auto-label variants `GraphName` for `graph_id` and `GraphGroupMasterName` for `graph_group_id`.
- Treat group-local orientation and LMB auto-resolution as dependent on singleton graph-group context.

I asked for clarification on whether `graph_group_master_id` was really intended, since group-master ids are not the same thing as graph-group ids.

## 2026-04-04 User

Do another deep review of the plan and re-reading all relevant part of the code to identify if there is any remaining point requiring clarification

## 2026-04-04 Codex

Deep-review planning response:

- The remaining design was tightened so `graph_group_id` is always populated whenever `graph_id` is known.
- Singleton `graph_id` selectors can provide graph-group context for group-local auto-label and auto-range resolution through the static graph-to-group mapping.
- No further product-level clarification appeared necessary from the current code shape.

## 2026-04-04 User

The integral from this card no longer converge to the same right value from before your changes.
But also:

`./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate generate integrate -c "quit -o"`

It crashes with

`selector = "value_range" is required for a value-range selector for key "default.selectors.graph_0_only" in TOML source string`

Why would a value range be necessary there?

## 2026-04-04 Codex

Investigation and fix summary:

- The `value_range` error was a flat-selector deserialization bug, not an intended requirement.
- The root cause was that the flat untagged selector helper could deserialize `selector = "discrete_range"` into the value-range helper first because their shapes overlapped.
- I fixed this by splitting the helper tag into selector-specific tags so `discrete_range` now deserializes correctly.
- I added a regression test for flat discrete-range selector deserialization and verified it with focused `cargo test` and `cargo check`.

## 2026-04-04 User

Add a graph_id, orientation_id (for each of the two graphs) and lmb channels (for each of the two graphs) discrete histograms to the epem_a_httx/LO example toml card, using all of the advanced setup we discussed, and also adjust the sampling settings in that example card so as to sample with MonteCarlo the two graphs, the orientations and the loop momentum channels.

## 2026-04-04 Codex

Implementation summary:

- Updated the `epem_a_ttxh/LO` example card to Monte Carlo sample graphs, orientations, and LMB channels.
- Added LO quantities `graph_id`, `orientation_id`, and `lmb_channel_id`.
- Added inactive singleton selectors `graph_0_only` and `graph_1_only`.
- Added discrete observables:
  - `graph_id` with auto graph-id domain and graph-name labels
  - `graph_0_orientation_id` / `graph_1_orientation_id`
  - `graph_0_lmb_channel_id` / `graph_1_lmb_channel_id`

## 2026-04-04 User

The integral from this card no longer converge to the same right value from before your changes.
But also:

`./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate -c "inspect -p epem_a_tth -i LO -x 0.1 0.2 0.3 0.4 0.5 0.6; quit -o"`

does not generate events anymore.
So you should investigate, and in particular check if the selectors still work as expected.

## 2026-04-04 Codex

Investigation summary:

- `inspect` not showing events was a separate CLI issue: that command path was not forcing `generate_events = true`. This is now fixed.
- The new named inactive selectors behave correctly:
  - they are ignored by runtime event rejection
  - they still gate the corresponding observables
- With explicit discrete Monte Carlo sampling active, the event metadata is populated correctly for `graph_id`, `graph_group_id`, `orientation_id`, and `lmb_channel_id`, and the per-graph orientation/LMB observables fill only on the matching graph selection.
- For the same seed and sample count, integrating with the new observables/selectors gives the same result as integrating without them.
- The remaining difference relative to the pre-change example comes from the changed sampling mode itself:
  - summed sampling at `20k` samples gave `-1.59(10)e-3` on the inspected phase
  - Monte Carlo sampling over graph/orientation/LMB at `20k` samples gave `-1.63(24)e-3`
  - these are statistically consistent, but the Monte Carlo discrete sampling has significantly larger variance, so convergence is slower/noisier in this example

The resulting plan is tracked in `discrete_histograms.md`.

## 2026-04-04 User

PLEASE IMPLEMENT THIS PLAN:

- one shared numeric `CutInfo` for runtime and output
- tagged continuous / discrete histogram settings and snapshots
- `graph_id`, `graph_group_id`, `orientation_id`, and `lmb_channel_id` as discrete quantities
- named selector references from observables
- selectors with `active: bool` plus an inclusive `discrete_range` variant
- auto domains and auto labels resolved once at runtime initialization

## 2026-04-04 Codex

Implementation start:

- Synced the tracking documents to the final approved plan.
- Starting with the selector/settings refactor because it drives quantity compilation, observable wiring, CLI/schema updates, and the later histogram-core extraction.

## 2026-04-04 Codex

Implementation status update:

- Implemented the tagged continuous/discrete histogram stack, including discrete bin ids, labels, ordering, merge semantics, snapshotting, JSON, and HwU output.
- Added numeric event metadata fields and discrete metadata-backed quantities for `graph_id`, `graph_group_id`, `orientation_id`, and `lmb_channel_id`.
- Added selector `active`, selector `discrete_range`, and observable named `selections`, with observable-local selector reuse.
- Exposed standalone histogram accumulation through the pyo3 API and updated Python snapshot/event metadata types.
- Updated schema/completion/template handling so selectors and observables keep their flat card shape while supporting the new discrete model.
- Migrated checked-in `cross_section` and `jet_count` observable cards to explicit discrete histogram definitions.
- Verified the implementation with formatting, compile checks, focused tests, and clippy.

## 2026-04-04 User

Add a `graph_id`, `orientation_id` (for each of the two graphs) and `lmb_channel_id` (for each of the two graphs) discrete histograms to the `epem_a_ttxh/LO` example TOML card, using the advanced selector-based setup, and adjust the sampling settings in that example so graphs, orientations, and loop-momentum channels are all Monte Carlo sampled.

## 2026-04-04 Codex

Follow-up example-card update:

- Updated `examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml`.
- Switched both `LO` and `LO_semi_inclusive` sampling blocks to Monte Carlo over:
  - graphs
  - orientations
  - LMB channels
- Added LO quantity definitions for:
  - `graph_id`
  - `orientation_id`
  - `lmb_channel_id`
- Added inactive named selectors:
  - `graph_0_only`
  - `graph_1_only`
  using singleton `graph_id` discrete ranges.
- Added advanced discrete observables:
  - one `graph_id` histogram with `GraphIds` + `GraphName`
  - one `orientation_id` histogram per graph with `OrientationIds` + `Orientation`
  - one `lmb_channel_id` histogram per graph with `LmbChannelIds` + `LmbChannelEdgeIds`
- Verified the example card remains valid TOML after the edit.

## 2026-04-04 User

Reported runtime failure when parsing an inactive `graph_id` selector with `selector = "discrete_range"`, receiving the error:

`selector = "value_range" is required for a value-range selector for key "default.selectors.graph_0_only" in TOML source string`

## 2026-04-04 Codex

Bugfix:

- Identified that the flat selector deserializer was using an untagged helper enum where `discrete_range` still structurally matched the `value_range` helper before the manual tag check ran.
- Tightened the flat selector serde helpers so each helper struct accepts only its own selector tag.
- Added an explicit deserialization test for flat `discrete_range` selectors.

## 2026-04-04 User

Reported that the `epem_a_ttxh/LO` example no longer integrated to the published targets and asked for deep debugging, with an explicit request to:

- first commit and push the current discrete-histogram work on a new `discrete_histograms` branch
- then verify whether the issue already existed before the discrete-histogram changes
- and if not, debug the new changes

## 2026-04-04 Codex

Saved and started the regression investigation:

- committed the discrete-histogram work locally with a dedicated commit
- pushed `discrete_histograms` to the remote and set upstream tracking
- reproduced the reported near-zero real-phase LO result on the current branch
- verified that the same issue already exists on the parent revision immediately before the discrete-histogram commit
- established that the physical LO value was living in the imaginary component before the discrete-histogram work

## 2026-04-04 User

Requested further investigation after confirming that the exact command

- `./gammaloop --clean-state ./examples/cli/epem_a_ttxh/LO/epem_a_tth_LO.toml run generate generate integrate -c "quit -o"`

returned near-zero real results instead of the published LO targets.

## 2026-04-04 Codex

Root cause and fix:

- narrowed the problem to the LO cross-section bare-cut assembly, not to selectors, observables, histogram filling, or the integration accumulator
- identified that the bare contribution built from `pass_two_result` was missing the subset-dependent complex phase
- restored that phase in
  - `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`
- rebuilt the CLI and verified:
  - with a reduced `20k`-sample run and `--show-phase both`, the LO result moved from the imaginary component to the real component
- with the exact user command, the real-phase results became
  - `epem_a_tth@LO = +1.547(21)e-3`
  - `epem_a_tth@LO_semi_inclusive = +5.330(72)e-4`
  - which are both within about `1σ` of the card targets

## 2026-04-04 User

Reported that with

```toml
[default_runtime_settings.integrator.observables_output]
format = "hwu"
```

the final emitted-results summary still listed only the JSON observable file, and requested that:

- the observable output mode support a list such as `["hwu", "json"]`
- a scalar like `"hwu"` become exclusive, not implicitly `"json"` as well
- the end-of-run emitted-files summary list all observable files when several are emitted

## 2026-04-04 Codex

Implemented the observable-output-format update:

- changed `ObservablesOutputSettings.format` from a single format to a list of formats
- initially accepted both scalar and list TOML forms
- updated the integration summary plumbing so each emitted observable file gets its own row in the final "Integration results emitted" table
- added tests covering:
  - scalar-or-list parsing
  - list parsing and round-tripping
  - summary rendering of multiple emitted observable files
- rebuilt the CLI and verified end to end:
  - with a temporary copy of the LO card using `format = ["hwu"]`, the final summary listed only
    - `integrands/epem_a_tth@LO/observables_final.hwu`
  - with a temporary copy of that card using `format = ["hwu", "json"]`, the final summary lists both
    - `integrands/epem_a_tth@LO/observables_final.hwu`
    - `integrands/epem_a_tth@LO/observables_final.json`

## 2026-04-04 User

Requested a follow-up simplification:

- remove support for the scalar form entirely
- only support list syntax, even for a single output format
- remove tests that try to verify support for a string value

## 2026-04-04 Codex

Applied the list-only follow-up:

- removed the custom scalar-or-list serde/schema compatibility layer from `ObservablesOutputSettings`
- the setting now accepts only list syntax such as:
  - `format = ["hwu"]`
  - `format = ["hwu", "json"]`
- replaced the old scalar-support unit test with a single-entry-list round-trip test
- migrated checked-in TOML files and example cards that still used scalar observable-output formats to list syntax
- updated the differential integration tests to pass list syntax through `set process kv`
- verified:
  - `cargo fmt`
  - `cargo test -p gammalooprs test_observables_output_settings_serialize_deserialize -- --nocapture`
  - `cargo test -p gammalooprs observables_output_settings_accept_single_entry_format_lists -- --nocapture`
  - `cargo test -p gammalooprs observables_output_settings_accept_multiple_formats -- --nocapture`
  - `cargo test -p gammalooprs results_output_summary_lists_all_emitted_observable_files -- --nocapture`
  - `cargo check -p gammaloop-api --features ufo_support`
- attempted to run the broader `gammaloop-integration-tests` differential target as well, but that target currently fails earlier due to an unrelated pre-existing compile error in `tests/tests/test_runs/events.rs` from a missing `force_generate_events` field in an `EvaluateSamples` initializer

## 2026-04-04 User

Requested three follow-up changes:

- verify that `--read-only-state` never lets the CLI write into the active state folder, and make such attempts fail cleanly
- in `--read-only-state`, move the default integration workspace out of the state and into the current working directory as `integration_workspace_<STATENAME>` or `integration_workspace`
- fix `run.toml` so it records all actually executed commands from boot-card commands, `run ...` command blocks, and inline `-c "..."`
- add Python API getters on a loaded session:
  - `get_run_history()`
  - `get_global_settings()`
  - `get_active_command_blocks()`

## 2026-04-04 Codex

Implemented the read-only-state and run-history follow-up:

- added a shared path-based write guard for the active state folder in `crates/gammaloop-api/src/lib.rs`
- applied it to:
  - `--clean-state`
  - integration workspace creation and updates
  - `save state`
  - `save dot`
  - `save standalone`
  - evaluator compilation triggered by `generate` when that would write under the active state
- changed the default read-only integration workspace location to:
  - `./integration_workspace_<state_name>`
- flattened run-history recording so the persisted `run.toml` captures the commands that actually executed rather than the wrapper `run ...` command itself
- added session-side string/dictionary getters and exposed them through the Python API as:
  - `get_run_history()`
  - `get_global_settings()`
  - `get_active_command_blocks()`
- verified with targeted Rust tests and the CLI integration test:
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
- initially the Python runtime validation was blocked because the editable extension in `.venv` was stale even though `cargo check -p gammaloop-api --features ufo_support,python_api` was already passing

## 2026-04-04 User

Commented that the Python rebuild should be achievable directly through the `just build-api` workflow, and that it appeared to have been figured out.

## 2026-04-04 Codex

Completed the Python runtime validation using the standard build path:

- cleared stale `maturin` / `cargo` build processes that were holding the target lock from earlier interrupted attempts
- reran `just build-api` successfully
- confirmed from `.venv` that the live module exposes:
  - `get_run_history`
  - `get_global_settings`
  - `get_active_command_blocks`
- ran the targeted Python API integration test successfully:
  - `PYTHON=/Users/vjhirsch/Documents/Work/gammaloop_main/.venv/bin/python cargo test -p gammaloop-integration-tests --features python-api-tests --test test_python_api python_session_getters_return_live_toml_and_command_blocks -- --nocapture`

## 2026-04-05 User

Requested a small CLI display adjustment:

- when the CLI prints `Running command <cmd>`
- if `<cmd>` is multiline and has more than 5 lines
- truncate the remainder and show `[...]` instead

## 2026-04-05 Codex

Implemented the command-display truncation in `crates/gammaloop-api/src/session.rs`:

- kept the change local to the shared `display_command(...)` path used by the `Running command` log lines
- multiline commands with more than 5 lines now display only the first 5 lines, followed by `[...]`
- single-line commands and multiline commands up to 5 lines are unchanged
- added unit tests covering both:
  - the exact 5-line boundary
  - truncation for commands with 6 lines
- verified with:
  - `cargo test -p gammaloop-api display_command_ -- --nocapture`
  - `cargo fmt`
  - `cargo check -p gammaloop-api`

## 2026-04-05 User

Requested a PR-review follow-up workflow:

- pull all comments from PR `#57` using `gh`
- address the relevant ones locally
- commit and push the fixes
- update the PR and reply briefly to each comment with the resolution

## 2026-04-05 Codex

Pulled the PR review threads with `gh` and addressed the actionable comments:

- fixed the discrete selector matching branch to match on a reference explicitly
- added explicit invalid-range checks on the Rust and Python discrete histogram constructors before allocation
- changed graph-group-context resolution to use the referenced quantity type instead of literal quantity names
- validated negative and out-of-range singleton graph/group selector values during graph-group-context resolution
- replaced unchecked label-resolution indexing with checked access and clear configuration errors
- added `DiscreteBinOrdering::as_str()` and used it in Python snapshot serialization so ordering strings stay in accepted snake_case form
- replaced the hot-path per-event graph-group scan with a cached graph-to-group mapping stored on the concrete process integrands
- cleaned the reviewed example-card nits in:
  - `examples/cli/gg_hhh/1L/gg_hhh_1L.toml`
  - `examples/cli/aa_aa/1L/aa_aa.toml`
- verified with:
  - `cargo fmt`
  - `cargo check`
  - `cargo test -p gammaloop-api display_command_ -- --nocapture`
  - `cargo test -p gammalooprs discrete_histogram_ -- --nocapture`
  - `cargo test -p gammalooprs resolve_graph_group_context_uses_selector_quantity_type -- --nocapture`
  - `cargo test -p gammalooprs resolve_graph_group_context_rejects_negative_singleton_values -- --nocapture`

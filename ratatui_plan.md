# Ratatui Integration Dashboard Plan

This file is the living plan for the ratatui integration-status dashboard work.
Update it whenever scope, decisions, or implementation sequencing changes.

## Confirmed decisions

- Ratatui becomes the default interactive streaming renderer via
  `--renderer=ratatui|tabled`, with `ratatui` as the new default.
- `info!()` output remains tabled/canonical. Ratatui is only for interactive
  streaming updates.
- If live updates are streamed but end-of-iteration summaries are not, the
  ratatui dashboard should suspend, let the tabled `info!()` summary print
  cleanly, and then resume on the next live update.
- The convergence chart uses:
  - x-axis = total sample count
  - y-axis = central value
- Multi-slot layout should use an all-slot ribbon plus a focused-slot detail
  area instead of a permanently wide all-columns matrix.
- ETA for the current iteration should be shown.
- Keyboard interrupt handling should become responsive inside
  `evaluate_samples(..)` / `evaluate_samples_raw(..)` after each individual
  sample finishes.

## Implementation status (2026-03-21)

### Implemented

- CLI renderer selection is live in
  `crates/gammaloop-api/src/commands/integrate.rs` with
  `--renderer=ratatui|tabled` and `ratatui` as the default.
- Ratatui streaming now runs through a dedicated dashboard controller with:
  - alternate-screen lifecycle
  - suspend/resume around canonical tabled `info!()` summaries
  - keyboard handling for tabs, slot focus, bin focus, sort mode, sort
    direction, metric toggles, chart-history window controls, chart-phase
    toggle, y-range controls/reset, help, `Ctrl-C`, and `x`
- Renderer-neutral status assembly was expanded with:
  - current-iteration timing metadata
  - training-slot / training-phase metadata
  - per-slot `chi^2` and max-weight impact
  - always-available statistics / max-weight / discrete sections for ratatui
- The ratatui dashboard currently includes:
  - overview tab with progress ribbon, ETA, convergence chart, focused-slot
    panel, results-summary table, and statistics composition bars
  - discrete tab with sortable/selectable bin rows and selected-bin detail
  - max-weight tab with overall and per-bin panels
- Ratatui now renders upstream `StyledText` content directly in tables and
  panels instead of flattening values to plain strings, so backend-independent
  formatting is shared across tabled and ratatui.
- Shared status-formatting helpers now cover ratatui-specific headers as well:
  - summary/discrete/detail headers for main results
  - slot-specific target-delta formatting
  - max-weight coordinate-table headers
  - statistics table rows and percentage-mix segments
- Shared value formatting now colors the parenthesized Monte Carlo uncertainty
  with the same threshold-based color as the corresponding relative error, and
  shared contribution headers now use the two-line `Contribution` /
  `(idx=...)` format.
- The canonical tabled renderer now explicitly flattens those contribution
  headers back to a single line, so tabled keeps
  `Contribution (idx=...)` while ratatui keeps the richer two-line layout.
- Shared component formatting now uses:
  - `re` in pink
  - `im` in yellow
  across both tabled and ratatui, including max-weight sign labels.
- Review-driven follow-ups now in place:
  - tabled raw-mode shutdown acknowledges only after terminal cleanup
  - Python API integrations stay on `tabled` until renderer selection is
    exposed there explicitly
  - the convergence history preserves the completed-iteration point when the
    next iteration emits its initial live update at the same sample count
  - `x` now discards a partial iteration instead of letting it affect later
    resume checkpoints
  - discrete-bin index sorting now respects ascending vs descending order
- The overview tab now uses the requested shape:
  - progress block with leading/trailing spacer rows, the reordered
    `[elapsed] | # samples total ... | Iteration ... | ETA ...` header, a
    right-aligned `#samples/s | /sample/core` group, and the gauge label kept
    on the bar
  - `Integrands` table with aligned `re:` / `im:` labels, focused-integrand
    marker, and optional per-integrand target displays using trimmed
    user-like scientific precision rather than raw noisy `f64` tails
  - large convergence panel with focused-integrand tracking, denser axis ticks,
    `current value +/- 4 sigma` y-window, chart-history controls that now
    scope by the actual iteration number seen in live updates instead of
    relying on completed-iteration summary markers, an interactive real/imag
    phase toggle, sigma-based y-range controls up to `±128σ`, and a title that
    shows the training-selection state plus the active history/y-range scope
  - the full-history view now preserves the global x-span by compacting old
    history points instead of simply dropping the earliest samples once the
    dashboard history buffer fills up
  - completed-iteration markers in the convergence plot are rendered as a
    distinct pink `iter` series in the legend so they remain visually separate
    from the live central-value line without overpowering it
  - focused-integrand table with persistent `% err` / `chi^2` / `m.w.i`
  - results-summary table with a dedicated second header row for per-integrand
    metric labels instead of repeating metric labels in each data row
- The statistics side now keeps the raw timing/eval/event rows readable in a
  tabled-like table and replaces the old ad hoc text blocks with proper mix
  panels for:
  - timing composition, sorted by largest slice first and split into
    `evaluators`, `itg_core`, `obs`, `param`, and `integrator`
  - precision mix, including the combined unstable/nan slice
- The discrete tab now follows the focused integrand for bin ordering/rendering
  and always shows `% err`, `chi^2`, and `m.w.i` in the discrete-bin table,
  plus `sample %`, while the selected-bin detail remains an all-integrands
  view.
- The discrete-tab row selection bug caused by the spacer row has been fixed so
  the highlighted bin and the selected-bin detail stay in sync.
- The max-weight tab now uses a wider coordinates panel filtered to the
  currently focused integrand, multi-line coordinate cells with explicit row
  heights, and the per-bin tables now follow the same contribution ordering as
  the discrete tab sorting for the first integrand. The top max-weight details
  table now also reserves enough height to keep the requested single trailing
  spacer row visible.
- Panel titles that mention a focused integrand now render the brackets in
  green and the integrand identifier in blue, to match the rest of the shared
  renderer formatting.
- The tabled renderer still obeys the existing show/hide flags; ratatui ignores
  those flags and always shows all available tabs/sections.
- Per-sample interrupt checks are implemented in:
  - `crates/gammalooprs/src/integrands/mod.rs`
  - `crates/gammalooprs/src/integrands/process/mod.rs`
- Ratatui now distinguishes:
  - `Ctrl-C` -> request full integration interrupt
  - `x` -> abort the current iteration after the current batch finishes
- Interactive streaming now shows a zero-progress frame immediately when each
  iteration starts, and the initial live status is now emitted before the
  expensive per-iteration worker-state construction begins.
- `--renderer=tabled` now starts a lightweight raw-mode control listener so
  streamed tabled output remains responsive to `Ctrl-C` and `x`.
- Tabled parity fixes now explicitly cover:
  - single-line contribution headers
  - preserved canonical `info!()` summaries while ratatui keeps its richer
    panel/table layout
- Targeted ratatui buffer tests were added in
  `crates/gammalooprs/src/integrate/mod.rs` for the overview, convergence-phase
  toggle, and discrete tabs.

### Still open / follow-up

- Python-facing renderer selection has not been threaded through every possible
  API surface yet; the CLI path is implemented.
- Interrupted work is still treated as partial-iteration work: the system
  remains checkpointed at completed-iteration boundaries only, which matches
  the current workspace/resume model. There is no partial-iteration persistence
  plan.
- Ratatui still needs broader state-machine coverage and more screenshot-driven
  polish:
  - mixed live ratatui + tabled summaries
  - non-TTY fallback
  - narrow-terminal layout behavior
  - max-weight-heavy payloads / long coordinate wrapping
- Optional UX polish still available:
  - more compact overview summary columns on narrower terminals
  - tighter chart/summary spacing on medium-width terminals
  - more aggressive chart-history downsampling tuning

## High-level goals

- Build a high-quality ratatui dashboard that feels native to terminal UIs
  rather than a direct port of the current tabled output.
- Preserve the existing cell content and label formatting conventions already
  captured in `crates/gammalooprs/src/integrate/display.rs`, and prefer using
  the preformatted display strings already carried by `StatusUpdate`.
- Keep the renderer boundary clean:
  - `gammalooprs` owns semantic integration status data
  - `gammaloop-api` owns TTY handling, alternate-screen lifecycle, and the
    ratatui event loop
- Improve readability for narrow terminals, multi-slot runs, and discrete-grid
  monitoring instead of forcing one static layout.

## Proposed architecture

### 1. Renderer selection and fallback behavior

- Add a `RendererOption` enum to
  `crates/gammaloop-api/src/commands/integrate.rs`.
- Thread that selection through CLI and Python bindings so both frontends stay
  aligned.
- Behavior matrix:
  - `--renderer=tabled`: keep the current streaming-in-place text renderer,
    plus keyboard controls for `Ctrl-C` and `x`.
  - `--renderer=ratatui` on a TTY: use the dashboard.
  - `--renderer=ratatui` on a non-TTY: log one fallback message and use current
    tabled logging behavior.
  - `--show-summary-only`: bypass ratatui entirely and print the current final
    tabled summary.

### 2. Semantic status snapshot refactor

- Keep `StatusUpdate` as the renderer-neutral contract, but make it less
  table-shaped and more complete for ratatui.
- Extend the status payload with data that ratatui needs but tabled currently
  hides or compresses:
  - explicit per-slot metrics for all slots, not just slot 0
  - current iteration timing metadata needed for ETA
  - training-slot and training-phase metadata for the live chart
  - typed statistics/timing breakdowns instead of only a prebuilt statistics
    table
  - monitored discrete-axis/path metadata for breadcrumbs and per-bin drilldown
- Prefer attaching behavior to existing section/state structs with methods where
  that improves discoverability, instead of adding more free helper functions.

### 3. Per-slot metrics model

- The current main table already carries per-slot value and relative error, but
  `chi^2` and max-weight impact are effectively slot-0 metadata.
- Ratatui should have access to per-slot:
  - value +/- uncertainty
  - relative error
  - `chi^2/dof` where meaningful
  - max-weight impact
  - target deltas for the first slot when targets exist
- Tabled can continue to render the compact slot-0 metadata columns, while the
  ratatui overview can show a richer all-slot metrics matrix.

### 4. Ratatui controller/event-loop design

- Implement ratatui in `gammaloop-api`, likely as a small dashboard module
  split into:
  - app state
  - terminal lifecycle/suspend-resume guard
  - drawing/layout helpers
  - key handling
  - status-history/cache logic
- Use a dedicated UI controller/thread:
  - integration sends `StatusUpdate` messages
  - UI thread handles input, redraw ticks, resize, suspend/resume, and shutdown
- This keeps the integrator hot path simple and makes interactive features like
  sorting and column toggles straightforward.

### 5. Suspend/resume lifecycle

- Ratatui path should use alternate screen + raw mode while active.
- Before emitting a tabled iteration summary or final summary:
  - leave alternate screen
  - disable raw mode
  - clear ratatui state from the terminal
  - print the tabled summary with the existing `info!()` path
- On the next live update, re-enter alternate screen and redraw from the latest
  cached dashboard state.

## Dashboard UX plan

### Global shell

- Keep a stable top-level frame with:
  - top title/tabs row
  - progress/status ribbon
  - main content area
  - bottom help/status footer
- Provide a small persistent key-hint footer, with a fuller help overlay on `?`.
- Support dynamic layout tiers:
  - wide: side-by-side panels
  - medium: stacked secondary panels
  - narrow: focus-first layout with fewer simultaneous tables

### Tab 1: Overview

- This is the main landing tab.
- Layout target on wide terminals:
  - top: progress ribbon + slot ribbon
  - left: convergence chart
  - right: focused-slot details card
  - bottom-left: all-slot metrics matrix
  - bottom-right: integration statistics panel

#### Progress ribbon

- Fancy ratatui gauge for current iteration progress.
- Include:
  - iteration number/state
  - completed/target iteration samples
  - total samples
  - throughput
  - ETA for the current iteration
  - elapsed wall time
- Use color states for healthy / warning / interrupted conditions.

#### Slot ribbon

- Show every slot as a compact card so multi-slot runs stay visible together.
- Each card should include at least:
  - slot label
  - current value +/- uncertainty
  - relative error badge
  - `chi^2/dof` badge
  - max-weight impact badge
- One slot is focused at a time; focus drives the right-side detail card and,
  when appropriate, secondary annotations in the chart.

#### Convergence chart

- Plot the training-slot integral in the selected training phase.
- Datasets:
  - central-value line
  - upper uncertainty line
  - lower uncertainty line
  - completed-iteration markers
  - target line when available
- Keep completed-iteration points visually distinct from intra-iteration live
  points.
- Use axis hysteresis/padding so the plot does not jitter on every redraw.
- Downsample history to terminal width.
- Allow:
  - phase toggle between real/imag for the focused integrand
  - history window toggle between full history and the last `N` iterations
  - y-range adjustment/reset in units of the current MC uncertainty

#### All-slot metrics matrix

- Show all slots together in a compact comparison matrix.
- Include interactively toggled extra columns:
  - relative error
  - `chi^2/dof`
  - max-weight impact
- Proposed interaction:
  - `r` toggles relative-error column
  - `c` toggles `chi^2/dof`
  - `w` toggles max-weight impact
  - `m` cycles density presets (`compact`, `metrics`, `full`)
- On narrow terminals, collapse this matrix into per-slot stacked rows instead
  of truncating columns aggressively.

#### Statistics panel

- Keep raw numbers readable, but complement them with visual composition bars.
- Show both exact values and percent-based visuals for metrics that should sum
  to 100%.
- Best candidates:
  - timing breakdown:
    - parameterization
    - evaluator
    - integrand
    - event processing
    - integrator overhead
  - precision mix:
    - double
    - quad
    - arb
  - stability mix:
    - stable
    - unstable
    - nan
  - event/selection efficiency:
    - accepted
    - rejected
- Preferred visual treatment:
  - compact numeric table on one side
  - horizontal composition bars on the other side
- That keeps the exact data visible while letting the user perceive balance and
  anomalies immediately.

### Tab 2: Discrete breakdown

- Always present, but show a clear empty-state panel when there is no monitored
  discrete dimension.
- Layout target:
  - left: discrete-bin table
  - right: selected-bin detail pane
  - optional bottom strip: slot comparison for the selected bin

#### Discrete-bin table

- Ratatui should improve on the current static sort choice with live sorting.
- Support:
  - sort mode cycle (`index`, `integral`, `error`, maybe `sample fraction`)
  - sort direction toggle
  - keyboard selection of a specific bin
- Show compact bar-style augmentations for:
  - sample fraction
  - target PDF
  - maybe absolute contribution magnitude

#### Selected-bin detail pane

- Show richer detail for the highlighted bin:
  - breadcrumb of fixed coordinates / monitored path
  - label formatting preserved from existing display helpers
  - per-slot values and uncertainties
  - per-slot relative errors
  - per-slot `chi^2/dof`
  - per-slot max-weight impact
  - max-weight coordinates if available

### Tab 3: Max weight

- Always present, with empty-state handling if no data exists.
- Layout target:
  - upper panel: overall max-weight details table
  - lower/side panel: per-bin max-weight details for either all bins or the
    currently selected discrete bin
- If a discrete tab selection exists, reuse that bin selection here so the UI
  feels coherent across tabs.
- For long coordinate payloads, show a wrapped viewer panel rather than forcing
  a huge wide table.

### Help/footer UX

- Bottom footer should show the current mode and the most important keys:
  - tab switching
  - slot/bin focus
  - sorting
  - metric toggles
  - density toggle
  - help
- `?` opens a centered help overlay/modal summarizing all keybindings.

## Creative ratatui features to leverage

- Composition bars for timing/precision/stability percentages.
- Focused-slot cards with threshold-colored metric badges.
- Empty-state panels that explain why a tab has no data instead of showing a
  blank area.
- Density modes so the same dashboard can feel spacious on large terminals and
  still usable on laptops.
- Mode/status toast line in the footer when a sort mode or column toggle
  changes.
- Stable selection coupling between tabs, especially for the focused slot and
  selected discrete bin.

## ETA plan

- Add explicit current-iteration elapsed time to the status snapshot so ETA does
  not need to be guessed from total integration runtime.
- Compute ETA from smoothed current-iteration throughput rather than using only
  the latest batch delta.
- Show:
  - instantaneous-ish throughput
  - smoothed throughput
  - ETA to iteration completion
- If too little data exists, show `ETA: warming up` instead of a noisy guess.

## Interrupt-responsiveness plan

- Add interrupt checks after each completed sample inside:
  - `crates/gammalooprs/src/integrands/mod.rs`
  - `crates/gammalooprs/src/integrands/process/mod.rs`
- Also break before starting the next sample once the interrupt flag is set, so
  no extra sample is started after the user has already interrupted.
- Do not change the completed-iteration checkpoint model: interrupted work
  inside the current iteration is still not persisted as a resumable partial
  iteration.
- Keep the existing outer batch-level checks in
  `crates/gammalooprs/src/integrate/mod.rs` as a second layer.
- Tests should verify prompt return at sample boundaries; persistence semantics
  remain iteration-level, not partial-iteration-level.

## Proposed implementation phases

### Phase 1: Renderer plumbing and status refactor

- Status: mostly done for the CLI/runtime path.
- Add renderer selection to CLI/Python/config surfaces.
- Add ratatui dependency to `crates/gammaloop-api/Cargo.toml`.
- Refactor `StatusUpdate` into a richer renderer-neutral snapshot.
- Add typed statistics/timing breakdown structures.
- Add per-slot metric support for all slots.

### Phase 2: Dashboard controller and terminal lifecycle

- Status: done for the current CLI dashboard implementation.
- Implement alternate-screen/raw-mode guard.
- Implement UI thread/controller with message passing.
- Add suspend/resume around iteration-summary and final-summary tabled output.
- Add resize handling and redraw ticks.

### Phase 3: Overview tab

- Status: done for the first production cut.
- Implement progress ribbon with ETA.
- Implement slot ribbon and focused-slot detail card.
- Implement convergence chart and status history cache.
- Add convergence controls for phase toggle, history-window selection, and
  y-range adjustment/reset in units of the current MC error.
- Implement all-slot metrics matrix with interactive metric toggles.
- Implement statistics table + composition bars.

### Phase 4: Discrete and max-weight tabs

- Status: done for the first production cut.
- Implement discrete-bin table with keyboard sorting and selection.
- Implement selected-bin detail pane.
- Implement overall and per-bin max-weight panels.
- Couple slot/bin focus across tabs where it helps.

### Phase 5: Interrupt responsiveness

- Status: core interrupt checks are done; iteration-level persistence semantics
  stay unchanged.
- Add per-sample interrupt checks in `evaluate_samples_raw` paths.
- Ensure partial batch results are returned cleanly.
- Add targeted tests.

### Phase 6: Polish and tests

- Status: in progress.
- Ratatui buffer/snapshot tests for wide and narrow layouts.
- Tabled parity tests to ensure the refactor does not regress the current
  non-ratatui summaries, including preserving single-line contribution headers
  in the canonical tabled output while ratatui keeps the richer two-line
  header layout.
- State-machine tests for:
  - live-only streaming
  - iteration-summary-only streaming
  - mixed ratatui live + tabled iteration summaries
  - final-summary suspend/teardown
  - non-TTY fallback
- Manual validation passes for:
  - single-slot
  - multi-slot
  - discrete monitoring
  - long max-weight coordinate payloads

### Phase 7: Accuracy-target termination and ETA-to-target

- Status: in progress.
- Add optional `integrator.target_relative_accuracy` and
  `integrator.target_absolute_accuracy` settings with `None` defaults and
  correct `SHOWDEFAULTS` serialization.
- Keep `integrator.n_max` as the hard stop, but also terminate after a
  completed iteration when either configured accuracy target is reached.
- Use only slot 0 and only slot 0's trained phase when evaluating accuracy
  targets, even in multi-slot runs.
- Add `ETA to target` to the ratatui overview header, based on the current
  slot-0 trained-phase error, current throughput, and `1/sqrt(N)` error
  scaling, and include the active target condition in the header as
  `ETA to target (<spec>) <eta>`, with relative targets shown as percentages
  in fixed notation above `1e-4%` and scientific notation at or below that
  threshold.
- If both accuracy targets are configured, compute the ETA using whichever
  target would be reached first; if the relative reference is zero, treat the
  relative target as inactive.

## Risks to watch

- Overfitting the semantic snapshot to the current tabled structure instead of
  giving ratatui the richer data it needs.
- Letting narrow-terminal behavior degrade into clipped unreadable tables.
- Recomputing expensive history or derived metrics on every frame instead of
  caching in the UI controller.
- Making suspend/resume flaky around `info!()` logging; this needs clean
  terminal ownership boundaries.
- Breaking current tabled output while refactoring shared status assembly.

#import "../products/shared.typ": source-link

= Integration dashboard architecture
<integration-dashboard-architecture>
#quote(block: true)[
#strong[Reviewed:] 2026-08-18 against the current command, event-loop, renderer, and renderer-test
sources

#strong[Lifecycle:] Current implementation contract. The earlier implementation plan and its
chronological status notes are retained separately as an archived record.
]

== User-facing boundary

The `integrate` command selects `ratatui` by default for live interactive updates and retains
`tabled` as the plain streaming renderer. `--no-stream-updates` disables in-iteration updates;
`--no-stream-iterations` disables end-of-iteration summaries. Both renderers consume the same
semantic `StatusUpdate` data, so choosing a terminal presentation must not change integration
state, convergence, results, or persisted artifacts.

The authoritative CLI arguments and renderer controller live in
#source-link("crates/gammaloop-api/src/commands/integrate.rs", label: "the integrate command").
Use the generated CLI reference for the exact flags accepted by a particular build.

== Ownership and data flow

The application layer owns terminal and command concerns:

+ `StreamingDisplayController` decides whether updates and iteration summaries are streamed.
+ `RatatuiStreamRenderer` sends semantic updates and lifecycle commands to one dashboard thread.
+ `RatatuiTerminal` owns raw mode, the alternate screen, cursor visibility, drawing, and cleanup.
+ `dashboard_event_loop` handles terminal input, redraws, suspend acknowledgements, and shutdown.

The domain layer owns renderer-independent integration status and the state needed to present
it. `RatatuiDashboardState` consumes `StatusUpdate`, retains bounded chart history, and renders
the overview, discrete-contribution, and max-weight views. Its implementation is in
#source-link("crates/gammalooprs/src/integrate/render_ratatui.rs", label: "the Ratatui renderer"),
while status assembly and renderer fixtures live in
#source-link("crates/gammalooprs/src/integrate/mod.rs", label: "the integration module").

```text
integration workers
  -> renderer-neutral StatusUpdate
  -> StreamingDisplayController
       -> tabled stream
       -> dashboard command channel
            -> terminal event loop
            -> RatatuiDashboardState
            -> stderr alternate screen
```

This boundary keeps terminal dependencies and input handling out of the numerical workers. The
dashboard is a view over status snapshots, not a second integration state machine.

== Terminal lifecycle and failure containment

The alternate screen is entered lazily on the first update. Before a canonical tabled summary
is printed, the controller synchronously suspends the dashboard and waits until the dashboard
thread has restored the terminal. Shutdown follows the same acknowledged path. Drop handlers
attempt cleanup as a final safeguard, and the event loop disables raw mode on exit.

`Ctrl-C` requests integration interruption; `x` requests aborting the current iteration. Other
keys only mutate dashboard view state: tab and slot selection, discrete-row sorting, metric
visibility, chart component/history/range, statistics scope, and help. View-state actions must
not mutate the scientific result.

== Verification contract

The maintained checks cover command parsing, keyboard action mapping, overview/slot metrics,
global versus focused statistics, chart phase and history behavior, ETA display, discrete-bin
detail, max-weight formatting, and narrow terminal layouts. Changes to status fields must update
both renderers and their fixtures; changes to terminal lifecycle must retain cleanup on normal
shutdown, suspension, interruption, and error paths.

The old plan proposed features and recorded intermediate states. It is available as the
#link("ratatui-integration-dashboard-history.typ")[archived dashboard implementation history],
but it is not an API or freshness authority.

#import "../../shared.typ": boundary, callout, developer-link, source-link

#let events = [
= Inspect events and observables

GammaLoop can evaluate selected points without launching a full adaptive integration. This is
the shortest path for inspecting the generated cut events, testing selectors, and checking that
histograms receive the expected weights. It is an investigative workflow: a few chosen points
cannot establish convergence or replace the persisted result of an integration.

== Start from the maintained differential example

The repository contains matching Python and Rust journeys for the differential Local Unitarity
process `e+ e- > d d~ g`. Their run cards generate the real two-graph process and configure:

- event generation and additional event weights;
- jet transverse-momentum, jet-count, and down-quark-energy quantities;
- a leading-jet selector; and
- continuous and discrete histogram observables.

The run card's `display_named_settings_examples` block lets you inspect each named object before
evaluation:

// docs-example: syntax
```text
display quantities -p #0
display quantities -p #0 leading_jet_pt
display selectors -p #0 leading_jet_pt_cut
display observables -p #0 leading_jet_pt_hist
```

Use these displays as a preflight check. A missing quantity name belongs to configuration; an
empty event list belongs to event generation or selection; an unexpected histogram belongs to
the quantity, phase, binning, or weights. Keeping those stages separate is more informative than
debugging the final number alone.

== Run the Python journey

Build the native package, then run the maintained script from the repository root. It creates a
fresh state beside the script, generates the configured process, evaluates one integration-space
point, evaluates a two-point batch, and evaluates one momentum-space point:

// docs-example: syntax
```sh
nix develop
uv venv .venv
source .venv/bin/activate
just build-api
python examples/api/python/epem_a_ddxg_xs_LO/inspect_events.py
```

The exact workflow is kept in
#source-link("examples/api/python/epem_a_ddxg_xs_LO/inspect_events.py", label: "the Python event-inspection script")
and its
#source-link("examples/api/python/epem_a_ddxg_xs_LO/run.toml", label: "differential run card").
The core API sequence is:

// docs-example: compile
```python
from pathlib import Path

import numpy as np
from gammaloop import GammaLoopAPI

example = Path("examples/api/python/epem_a_ddxg_xs_LO")
api = GammaLoopAPI(
    state_folder=example / "state",
    boot_commands_path=example / "run.toml",
    clean_state=True,
)
api.run("run display_named_settings_examples")

point = [0.17, 0.31, 0.53, 0.23, 0.41, 0.67]
single = api.evaluate_sample(point, return_events=True)

points = np.array([
    point,
    [0.11, 0.29, 0.47, 0.19, 0.37, 0.59],
], dtype=float)
batch = api.evaluate_samples(points, return_events=True)

print(single.integrand_result, single.generated_event_count)
for group in single.event_groups:
    for event in group.events:
        print(event.cut_info.graph_id, event.weight, event.outgoing_momenta)
print(batch.observables["leading_jet_pt_hist"].sample_count)
```

`return_events=True` temporarily requests returned events for that call; it does not permanently
rewrite the integrand setting. Every `EventGroup` contains correlated accepted events produced by
one graph group. Each event carries incoming/outgoing momenta, its primary complex weight, cut and
graph identifiers, and any configured additional weights.

#callout("Single and batch results have different aggregation boundaries", [
  `evaluate_sample` returns one `EvaluationResult`: its event groups and observable snapshot
  belong to that sample. `evaluate_samples` returns per-row `SampleEvaluationResult` objects and
  one `BatchEvaluationResult.observables` dictionary merged across the whole batch. Do not read a
  batch histogram as if it belonged to the last row, and do not add the same batch snapshot once
  per sample.
])

The histogram values are immutable snapshots. Inspect `bins`, underflow/overflow, raw sums of
weights and squared weights, entry counts, and the aggregate statistics before deriving a plot.
Keep raw statistics when combining results; rounded rendered values are not a merge format.

== Run the Rust journey

The Rust example loads the same kind of run card through `StateLoadOption`, executes the display
block through a CLI session, then calls the typed evaluation functions. It is a `rust-script`
program so its dependency revision is recorded with the example:

// docs-example: syntax
```sh
NO_SYMBOLICA_OEM_LICENSE=1 \
EXTRA_MACOS_LIBS_FOR_GNU_GCC=T \
SYMBOLICA_LICENSE="$SYMBOLICA_LICENSE" \
rust-script --debug examples/api/rust/epem_a_ddxg_xs_LO/inspect_events.rs
```

The documentation gate compiles the request boundary below against the workspace API. It keeps
the event-return override and every required field synchronized even though the licensed process
evaluation remains outside the pull-request tier:

// docs-example: compile gammaloop-event-inspection
```rust
use gammaloop_api::commands::evaluate_samples::EvaluateSamples;
use ndarray::arr2;

let points = arr2(&[[0.17, 0.31, 0.53, 0.23, 0.41, 0.67]]);
let request = EvaluateSamples {
    process_id: None,
    integrand_name: None,
    use_arb_prec: false,
    minimal_output: false,
    return_generated_events: Some(true),
    momentum_space: false,
    points: points.view(),
    integrator_weights: None,
    discrete_dims: None,
    graph_names: None,
    orientations: None,
};
assert_eq!(request.return_generated_events, Some(true));

let loop_momenta = arr2(&[[
    0.11, -0.07, 0.19, // k1 = (px, py, pz)
    -0.13, 0.05, 0.29, // k2 = (px, py, pz)
]]);
let momentum_request = EvaluateSamples {
    process_id: None,
    integrand_name: None,
    use_arb_prec: true,
    minimal_output: true,
    return_generated_events: None,
    momentum_space: true,
    points: loop_momenta.view(),
    integrator_weights: None,
    discrete_dims: None,
    graph_names: None,
    orientations: None,
};
assert_eq!(momentum_request.points.ncols(), 2 * 3);
```

Read
#source-link("examples/api/rust/epem_a_ddxg_xs_LO/inspect_events.rs", label: "the Rust event-inspection script")
with its
#source-link("examples/api/rust/epem_a_ddxg_xs_LO/run.toml", label: "run card"). The ordinary
`EvaluateSamples` path returns the cross-language `f64` boundary. Rust callers that must retain
the active numeric precision use `EvaluateSamplesPrecise`; arbitrary-precision internal
evaluation does not make the ordinary result fields arbitrary-precision.
The historical `use_arb_prec` name denotes the compatibility override described below: it selects
configured `f128`, falling back to Arb only when `f128` is unavailable.
The #link("reference/interfaces/#sample-evaluation-contract")[sample evaluation contract]
defines this flattened momentum layout and the precision behavior once for CLI, Rust, and Python.

#boundary("Verification and resource cost", [
  The pull-request gate parses the commands, compiles the Python source, and type-checks the Rust
  request literal, but does not run process generation or native evaluation. The maintained
  scripts are the real acceptance path and require a built API, the repository's scientific
  dependencies, and an appropriate Symbolica license. They create or replace only their
  example-local `state` directory because they opt into `clean_state=True`; remove that option
  before pointing at work you need to keep.
])

== From point inspection to integration output

Selected points answer structural questions: whether events are produced, how graph/cut metadata
is grouped, which selector accepts them, and which histogram bins are filled. A production result
still needs an integration workspace, target-accuracy settings, iteration diagnostics, and the
final observable output. Continue with the #link("guides/process-generation/")[state and
integration-workspace guide] for resume semantics and the distinction between
`integration_result.json`, `observables_final.*`, and the authoritative checkpoint.

Contributors changing progress reporting, cancellation, or terminal ownership should also read
the #developer-link(
  "integration-dashboard",
  "ratatui-integration-dashboard.typ",
  "integration dashboard architecture",
).

The exact Python return types and fields are linked from the
#link("reference/interfaces/")[interface guide] and generated
#link("reference/python/gammaloop-python/")[Python API].
]

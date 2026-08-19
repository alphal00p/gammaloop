#import "../../shared.typ": boundary

#let api = [
= CLI, Rust, and Python APIs

== Rust packages

The Rust API spans two primary packages.

- `gammalooprs` contains the physics implementation, data model, integrands, integration,
  observables, and lower-level algorithms.
- `gammaloop-api` provides the supported state-loading API, session/command integration,
  CLI assembly, and the PyO3 extension used by the Python distribution.

Load a state with `StateLoadOption::load()`. Its options select a state directory, boot card,
logging overrides, read-only behavior, and optional settings overrides. Once loaded, callers
can create a CLI-style session or use dedicated structured operations. A raw command string is
the compatibility fallback, not a replacement for a typed operation where one exists.

```rust
use std::path::PathBuf;
use gammaloop_api::{state::CommandHistory, StateLoadOption};

let mut loaded = StateLoadOption {
    state_folder: Some(PathBuf::from("./state")),
    ..StateLoadOption::default()
}.load()?;

let command = CommandHistory::from_raw_string("display settings global")?;
loaded.cli_session().execute_command(command)?;
```

Some Rust types are public so that GammaLoop's packages can work together, but are not intended
as stable integration points. Prefer the state-loading and structured-operation APIs described
here; use the #link("reference/rust/")[Rust orientation] to choose a crate, then use that
revision's Rustdoc for exact signatures, trait implementations, and lower-level types.

=== Rust reference map

This authored map separates the supported workflow from the much larger compiled surface.
Follow these exact Rustdoc paths instead of beginning in an arbitrary internal module:

- *Load and command a state:* begin with
  #link("reference/rust/gammaloop_api/struct.StateLoadOption.html")[`StateLoadOption`],
  retain the returned
  #link("reference/rust/gammaloop_api/struct.LoadedState.html")[`LoadedState`], and
  execute parsed
  #link("reference/rust/gammaloop_api/state/struct.CommandHistory.html")[`CommandHistory`]
  values through
  #link("reference/rust/gammaloop_api/session/struct.CliSession.html")[`CliSession`].
- *Replay and inspect:* use
  #link("reference/rust/gammaloop_api/state/struct.RunHistory.html")[`RunHistory`] for run
  cards and
  #link("reference/rust/gammaloop_api/integrand_info/struct.IntegrandInfo.html")[`IntegrandInfo`] for
  structured graph-group, orientation, loop-basis, cut, and threshold metadata.
- *Evaluate samples:* construct
  #link("reference/rust/gammaloop_api/commands/evaluate_samples/struct.EvaluateSamples.html")[`EvaluateSamples`]
  for the ordinary `f64` boundary. Use
  #link("reference/rust/gammaloop_api/commands/evaluate_samples/struct.EvaluateSamplesPrecise.html")[`EvaluateSamplesPrecise`]
  only when the caller must retain the active `f64`, `f128`, or arbitrary-precision result type.
- *Work below the facade:* the `gammalooprs` Rustdoc includes the
  #link("reference/rust/gammalooprs/integrands/trait.HasIntegrand.html")[integrand contract],
  #link("reference/rust/gammalooprs/settings/struct.RuntimeSettings.html")[runtime settings],
  precision-specific result records, and the mergeable
  #link("reference/rust/gammalooprs/observables/struct.HistogramSnapshot.html")[histogram snapshot]
  boundary. Treat public items not named by a maintained workflow as implementation-facing
  unless their crate documents a stronger compatibility promise.

== Python packaging

#boundary("A standalone GammaLoop distribution", [
  The Python distribution and import package are both named `gammaloop`. Its compiled backend
  is `gammaloop._gammaloop`, which applications normally access through the public package.
  This is separate from the `symbolica.community.*` modules used by Spenso, Idenso, and Vakint;
  installing GammaLoop does not install those community modules as independent packages.
])

The main Python entry point is
#link("reference/python/gammaloop-python/GammaLoopAPI/")[`GammaLoopAPI`]. Its
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-new-associatedfunction")[constructor]
loads or creates one stateful session:

// docs-example: compile
```python
from gammaloop import GammaLoopAPI

gl = GammaLoopAPI(
    state_folder="./state",
    boot_commands_path="./run.toml",
)
gl.run("display settings global")
```

The Python package requires Python 3.11 or newer. Run `just build-api` when building the bindings
from a source checkout.

=== Python reference map

The generated module reference covers every registered public export, but the objects form a few
connected workflows rather than forty unrelated entry points:

- *Session lifecycle:* construct
  #link("reference/python/gammaloop-python/GammaLoopAPI/")[`GammaLoopAPI`], then use
  #link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-run-method")[`run`], the history
  getters, and
  #link("reference/python/gammaloop-python/SettingsValue/")[`SettingsValue`] to automate
  the same state and commands as the CLI.
- *Point evaluation:* `evaluate_sample` returns
  #link("reference/python/gammaloop-python/EvaluationResult/")[`EvaluationResult`], whose
  `sample` is a
  #link("reference/python/gammaloop-python/SampleEvaluationResult/")[`SampleEvaluationResult`].
  `evaluate_samples` returns
  #link("reference/python/gammaloop-python/BatchEvaluationResult/")[`BatchEvaluationResult`]
  with the same sample records and one batch-level observable snapshot.
- *Generated-integrand structure:* start at
  #link("reference/python/gammaloop-python/IntegrandInfo/")[`IntegrandInfo`], then follow its
  graph groups into graph, orientation, loop-momentum-basis, cut, and threshold records. These
  objects describe compiled structure; they do not mutate it.
- *Events and observables:* evaluation records lead to
  #link("reference/python/gammaloop-python/EventGroup/")[`EventGroup`] and
  #link("reference/python/gammaloop-python/Event/")[`Event`]. For caller-owned aggregation,
  #link("reference/python/gammaloop-python/HistogramAccumulator/")[`HistogramAccumulator`]
  produces immutable
  #link("reference/python/gammaloop-python/HistogramSnapshot/")[`HistogramSnapshot`]
  records with raw statistics that remain mergeable and reconstructible.

#boundary("Full integrations still use the command interface", [
  The current Python module does not expose a supported structured `integrate()` method. Run an
  integration through `GammaLoopAPI.run("integrate ...")` or the CLI and consume its persisted
  workspace/results. Several integration-result record classes are registered by the native
  module for ongoing API work, but without a public producer they are not part of the curated
  Python boundary yet.
])

#boundary("Result ownership", [
  Objects returned by evaluation and integration are snapshots. Mutating a Python list or
  dictionary derived from them does not update the live GammaLoop session. The explicit mutable
  exception is `HistogramAccumulator`, whose methods update caller-owned aggregation state.
])

== Inspect and evaluate an existing state

The following workflow assumes that `./state` contains a generated integrand. Evaluation points
are integrand-specific: their length must match the selected integrand and any discrete dimensions.
Use the generated reference for the exact
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-evaluate-sample-method")[single-sample]
and
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-evaluate-samples-method")[batch]
contracts.

// docs-example: compile
```python
from collections.abc import Sequence

from gammaloop import GammaLoopAPI

api = GammaLoopAPI(state_folder="./state", read_only_state=True)
api.run("display processes")

print(api.get_run_history())
runtime = api.get_default_runtime_settings()
print(runtime.get("integrator.n_start"))


def evaluate_one(point: Sequence[float]) -> complex:
    result = api.evaluate_sample(
        point,
        process_id=0,
        minimal_output=True,
    )
    return result.integrand_result
```

#boundary("Verification tier: compile", [
  The documentation checks compile this example as Python syntax. They do not import or execute
  the native module: a runtime check additionally needs a built GammaLoop package, a provisioned
  backend and license, an existing generated state, and a point with the correct dimension.
])

The #link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-run-method")[`run` method]
uses the same command language as the CLI and can change the in-memory state, settings, and run
history. `read_only_state=True` protects files inside the active state directory; it does not make
the Python object immutable or automatically persist the session. Inspect the live session through
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-run-history-method")[`get_run_history`],
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-global-settings-method")[`get_global_settings`],
and
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-active-command-blocks-method")[`get_active_command_blocks`].

For structured inspection,
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-integrand-info-method")[`get_integrand_info`]
describes the selected backend and graph structure, while
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-integrand-settings-method")[`get_integrand_settings`]
and
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-get-default-runtime-settings-method")[`get_default_runtime_settings`]
return detached, read-only `SettingsValue` snapshots. Use `get(path)`, attribute access, indexing,
or `to_dict()` to read them; modifying derived Python values does not update the live session.

== Sample evaluation and precision

#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-evaluate-sample-method")[`evaluate_sample`]
returns one sample result and the observable bundle for that one-sample batch.
#link("reference/python/gammaloop-python/GammaLoopAPI/#exports-gammaloopapi-evaluate-samples-method")[`evaluate_samples`]
accepts a two-dimensional NumPy array and returns per-sample results plus a batch-global observable
bundle. Per-sample weights, discrete coordinates, graph names, and orientations must have the same
row count as the input batch when provided. Graph and orientation selection applies only to
momentum-space evaluation.

Both the ordinary Rust and Python endpoints expose numeric results through an `f64` contract, even
when `use_arb_prec=True` requests arbitrary-precision internal evaluation. Evaluation may warm the
integrand and update in-memory caches or observable snapshots, including in a read-only-state
session. Rust-only callers that must retain the active precision use `evaluate_sample_precise` and
`evaluate_samples_precise`.

For a complete, source-backed run card and scripts that expose event groups, cut metadata,
selectors, and merged histogram snapshots, follow the
#link("guides/events-and-observables/")[events and observables guide].

== Symbolic conversion helpers

The supported module-level helpers are
#link("reference/python/gammaloop-python/atom_to_canonical_string/")[`atom_to_canonical_string`]
for parsing and canonicalizing a Symbolica expression,
#link("reference/python/gammaloop-python/evaluate_graph_overall_factor/")[`evaluate_graph_overall_factor`]
for evaluating a graph's symbolic prefactor, and
#link("reference/python/gammaloop-python/to_dots/")[`to_dots`] for rewriting tensor
contractions into Idenso dot-product notation. Their generated entries document accepted strings,
return values, failure modes, and task-oriented examples.

== CLI and settings reference

=== Shell completion

`--completions` emits a static script from the same Clap command tree as `--help`. Load it for
the current shell or save it in the shell's normal completion directory:

// docs-example: syntax
```sh
# Bash for the current session
source <(./gammaloop --completions bash)

# Zsh for the current session
source <(./gammaloop --completions zsh)
```

Fish, Elvish, PowerShell, and Nushell are also supported values. For Fish, pipe the output to
`source`; for Nushell, save the generated module and then source it. Bash and Fish output also
registers the repository wrapper spelling `./gammaloop`.

#boundary("Static and state-aware completion", [
  The generated shell script covers executable subcommands, flags, enumerated values, and path
  hints. Completion inside the interactive GammaLoop prompt additionally knows the active
  state's processes, integrands, settings, model objects, graphs, cuts, and orientations. A
  static shell script cannot contain that changing session data.
])

Use `./gammaloop --help` and the generated #link("reference/cli/")[CLI and settings reference]
for command names, aliases, flags, positional arguments, defaults, possible values, and setting
paths. The tutorials explain how commands and settings combine into a persistent workflow. The
#link("guides/diagnostics/")[logging and diagnostics guide] covers startup precedence, selective tag
filters, sink routing, release-build limitations, and copyable investigation patterns.
]

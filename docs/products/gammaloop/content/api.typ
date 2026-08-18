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
here; start with the #link("reference/rust/")[curated Rust API] and use its Rustdoc links when
you need lower-level types or trait implementations.

== Python packaging

#boundary("A standalone GammaLoop distribution", [
  The Python distribution and import package are both named `gammaloop`. Its compiled backend
  is `gammaloop._gammaloop`, which applications normally access through the public package.
  This is separate from the `symbolica.community.*` modules used by Spenso, Idenso, and Vakint;
  installing GammaLoop does not install those community modules as independent packages.
])

The main Python entry point is
#link("reference/python/gammaloop-python/#exports-gammaloopapi")[`GammaLoopAPI`]. Its
#link("reference/python/gammaloop-python/#exports-gammaloopapi-new-associatedfunction")[constructor]
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

== Inspect and evaluate an existing state

The following workflow assumes that `./state` contains a generated integrand. Evaluation points
are integrand-specific: their length must match the selected integrand and any discrete dimensions.
Use the generated reference for the exact
#link("reference/python/gammaloop-python/#exports-gammaloopapi-evaluate-sample-method")[single-sample]
and
#link("reference/python/gammaloop-python/#exports-gammaloopapi-evaluate-samples-method")[batch]
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

The #link("reference/python/gammaloop-python/#exports-gammaloopapi-run-method")[`run` method]
uses the same command language as the CLI and can change the in-memory state, settings, and run
history. `read_only_state=True` protects files inside the active state directory; it does not make
the Python object immutable or automatically persist the session. Inspect the live session through
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-run-history-method")[`get_run_history`],
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-global-settings-method")[`get_global_settings`],
and
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-active-command-blocks-method")[`get_active_command_blocks`].

For structured inspection,
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-integrand-info-method")[`get_integrand_info`]
describes the selected backend and graph structure, while
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-integrand-settings-method")[`get_integrand_settings`]
and
#link("reference/python/gammaloop-python/#exports-gammaloopapi-get-default-runtime-settings-method")[`get_default_runtime_settings`]
return detached, read-only `SettingsValue` snapshots. Use `get(path)`, attribute access, indexing,
or `to_dict()` to read them; modifying derived Python values does not update the live session.

== Sample evaluation and precision

#link("reference/python/gammaloop-python/#exports-gammaloopapi-evaluate-sample-method")[`evaluate_sample`]
returns one sample result and the observable bundle for that one-sample batch.
#link("reference/python/gammaloop-python/#exports-gammaloopapi-evaluate-samples-method")[`evaluate_samples`]
accepts a two-dimensional NumPy array and returns per-sample results plus a batch-global observable
bundle. Per-sample weights, discrete coordinates, graph names, and orientations must have the same
row count as the input batch when provided. Graph and orientation selection applies only to
momentum-space evaluation.

Both the ordinary Rust and Python endpoints expose numeric results through an `f64` contract, even
when `use_arb_prec=True` requests arbitrary-precision internal evaluation. Evaluation may warm the
integrand and update in-memory caches or observable snapshots, including in a read-only-state
session. Rust-only callers that must retain the active precision use `evaluate_sample_precise` and
`evaluate_samples_precise`.

== Symbolic conversion helpers

The supported module-level helpers are
#link("reference/python/gammaloop-python/#exports-atom-to-canonical-string")[`atom_to_canonical_string`]
for parsing and canonicalizing a Symbolica expression,
#link("reference/python/gammaloop-python/#exports-evaluate-graph-overall-factor")[`evaluate_graph_overall_factor`]
for evaluating a graph's symbolic prefactor, and
#link("reference/python/gammaloop-python/#exports-to-dots")[`to_dots`] for rewriting tensor
contractions into Idenso dot-product notation. Their generated entries document accepted strings,
return values, failure modes, and task-oriented examples.

== CLI and settings reference

Use `./gammaloop --help` and the generated #link("reference/cli/")[CLI and settings reference]
for command names, aliases, flags, positional arguments, defaults, possible values, and setting
paths. The tutorials explain how commands and settings combine into a persistent workflow. The
#link("manual/diagnostics/")[logging and diagnostics guide] covers startup precedence, selective tag
filters, sink routing, release-build limitations, and copyable investigation patterns.
]

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
here; consult the full Rust API reference when you need lower-level types or trait
implementations.

== Python packaging

#boundary("A standalone GammaLoop distribution", [
  The Python distribution and import package are both named `gammaloop`. Its compiled backend
  is `gammaloop._gammaloop`, which applications normally access through the public package.
  This is separate from the `symbolica.community.*` modules used by Spenso, Idenso, and Vakint;
  installing GammaLoop does not install those community modules as independent packages.
])

The main Python entry point is `GammaLoopAPI`:

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

== Sample evaluation and precision

`evaluate_sample` returns one sample result and the observable bundle for that one-sample
batch. `evaluate_samples` returns per-sample results plus a batch-global observable bundle.
Both the ordinary Rust and Python endpoints expose numeric results through an `f64` contract,
even if arbitrary precision was used internally. Rust-only callers that must retain the active
precision use `evaluate_sample_precise` and `evaluate_samples_precise`.

== CLI and settings reference

Use `./gammaloop --help` and the CLI reference for command names, aliases, flags, positional
arguments, defaults, and possible values. The settings reference lists available paths and
defaults. The tutorials explain how commands and settings combine into a persistent workflow.
]

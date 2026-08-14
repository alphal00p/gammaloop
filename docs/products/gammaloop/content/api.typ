#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust, Python, and CLI reference boundaries

#catalog-contract(
  rust-scope: "gammaloop",
  python-scope: "gammaloop",
)

== Rust ownership

The Rust surface is split across two primary packages.

- `gammalooprs` contains the physics implementation, data model, integrands, integration,
  observables, and lower-level algorithms.
- `gammaloop-api` owns the supported state-loading boundary, session/command integration,
  CLI assembly, and the PyO3 extension used by the Python distribution.

The stable conceptual entry point is `StateLoadOption::load()`. Load-time options select a
state directory, boot card, logging overrides, read-only behavior, and optional settings
overrides. Once loaded, callers can create a CLI-style session or use dedicated structured
operations. A raw command string is the compatibility fallback, not a replacement for a
typed operation where one exists.

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

Generated Rust reference pages must distinguish public supported catalog entries from
implementation-public types. Conventional Rustdoc remains available for the complete compiled
surface and trait implementations.

== Python packaging

#boundary("A standalone GammaLoop distribution", [
  The repository root `pyproject.toml` builds the Python distribution named `gammaloop` with
  Maturin. Users import the public package as `gammaloop`; its compiled implementation module
  is `gammaloop._gammaloop`. This is separate from the `symbolica.community.*` modules used by
  Spenso, Idenso, and Vakint. Building or installing GammaLoop does not install those community
  modules as independent Python packages.
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

The Python package requires Python 3.11 or newer. `just build-api` is the local development
build. Distribution builds may enable the stable ABI feature; that choice changes binary
compatibility, not the public module name.

== Sample evaluation and precision

`evaluate_sample` returns one sample result and the observable bundle for that one-sample
batch. `evaluate_samples` returns per-sample results plus a batch-global observable bundle.
Both the ordinary Rust and Python endpoints expose numeric results through an `f64` contract,
even if arbitrary precision was used internally. Rust-only callers that must retain the active
precision use `evaluate_sample_precise` and `evaluate_samples_precise`.

The generated reference must present these as distinct contracts. It must not suggest that the
Python API returns arbitrary-precision numeric values.

== CLI and settings metadata

The executable command tree is generated from Clap and settings descriptions are generated
from the existing schema/default machinery. These generated catalogs own command names,
aliases, flags, positionals, defaults, possible values, and settings paths. Hand-written guides
own sequencing, state effects, and examples.

Source starting points are #source-link("crates/gammaloop-api/src", label: "gammaloop-api") and
#source-link("crates/gammalooprs/src", label: "gammalooprs").
]

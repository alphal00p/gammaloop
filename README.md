<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-dark.svg">
  <img src="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-light.svg" width="300">
</picture>
</div>


*Computation of differential cross-sections using Local Unitarity.*

<img src="https://nix-ci.com/badge/gh:alphal00p:gammaloop/main" alt="NixCI Badge">


See [www.alphaloop.ch](https://www.alphaloop.ch) for the broader project and related literature.
See the [wiki](https://wiki.alphaloop.ch/) for command syntax and usage notes.
Process-generation syntax is documented here: [gammaLoop/ProcessGeneration](https://wiki.alphaloop.ch/en/gammaLoop/ProcessGeneration).
Implemented architecture notes live in [docs/architecture/architecture-current.md](docs/architecture/architecture-current.md).

## Quick start

GammaLoop is driven from a stateful CLI. A run card can initialize a state, execute commands, and leave behind a reusable working directory that can be resumed later.

Build the CLI from source:

```bash
just build-cli-release
```

The built binaries live under `./target/<profile>/gammaloop`. The repository also provides a thin `./gammaloop` wrapper for local development, so from the project root directory you can run:

```bash
./gammaloop ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml
```

This example:
- creates the state in `./examples/cli/gg_hhh/1L/state`
- imports the Standard Model
- generates one pentagon contribution to the `gg -> hhh` one-loop amplitude
- integrates a kinematic configuration without thresholds
- saves the resulting state on `quit -o`

You can then bootstrap from that existing state and run additional command blocks:

```bash
./gammaloop -s ./examples/cli/gg_hhh/1L/state run integrate_physical -c "quit -o"
```

That command reuses the existing state, loads the command-block definitions stored in `./examples/cli/gg_hhh/1L/state/run.toml`, executes the `integrate_physical` block defined there, and persists the updated state on exit.

If you also want to render the generated diagrams, change to `./examples/cli/gg_hhh/1L/state` and run:

```bash
just draw
```

## Installation from source

GammaLoop is currently developed primarily from source.

Requirements:
- Rust `1.85+`
- `just` (a lightweight make-like command runner; install with `cargo install just`)
- a recent GNU toolchain
- Python `3.11+` only if you also want the Python bindings
- FORM `v4.2.1+` (only used for analytical integration of integrated UV counterterms)
- Python module [ufo-model-loader](https://pypi.org/project/ufo-model-loader/) if you want gammaloop to be able to import UFO models

Additional requirements for rendering diagrams:

```bash
cargo install clinnet
cargo install typst
```

If local dependency management is inconvenient, you can instead use the provided Nix flake:

```bash
nix develop
```
which will bring you to an environment with all dependencies available.

Common build commands:

```bash
just build-cli
just build-cli-release
just build-api
just fmt
just clippy
just test TEST_NAME
```

Notes:
- `just build-cli` builds the CLI in the `dev-optim` profile at `./target/dev-optim/gammaloop`
- `just build-cli-release` builds the release CLI at `./target/release/gammaloop`
- the repo-root `./gammaloop` wrapper executes the built binary for you during development

## Shell completion

The CLI can emit static shell-completion scripts from its Clap command definition that can be loaded into your shell as follows:

```bash
# bash
source <(./gammaloop --completions bash)
# zsh
source <(./gammaloop --completions zsh)
# fish
./gammaloop --completions fish | source
# nushell
./gammaloop --completions nushell | save --force /tmp/gammaloop-completions.nu
source /tmp/gammaloop-completions.nu
```

These completions cover the executable command surface (subcommands, flags, enum values, and path hints). The richer state-aware completion remains available inside the GammaLoop REPL. Bash and Fish also bind completion to the repository wrapper path `./gammaloop`.

## CLI model

The current CLI revolves around a persistent state folder.

- start from a run card: `gammaloop <card.toml>`
- resume an existing state: `gammaloop -s <STATE_FOLDER>`
- run command blocks from the active state: `gammaloop -s <STATE_FOLDER> run <block_name>`
- add ad-hoc commands to a run invocation with `run ... -c "..."`
- define new command blocks interactively with `start_commands_block` / `finish_commands_block`

By default the state contains, among other files:
- `global_settings.toml`
- `default_runtime_settings.toml`
- `run.toml` — a high-level audit log of the session, including the frozen boot settings, command-block definitions, and the top-level commands that built the current state
- `processes/`

In particular, `run.toml` is intended to be replayable: running

```bash
./gammaloop <STATE_FOLDER>/run.toml
```

should rerun the recorded commands and rebuild the current state from scratch to the same point.

## Rust and Python API

The CLI remains the primary interface, but the same state-loading entry point is also exposed from Rust and Python. In both cases the entry point only handles static load-time options such as:
- state folder
- boot card
- logging overrides
- read-only mode
- optional settings override files

Once a state is loaded, there are dedicated API methods available for some tasks. The generic fallback is to run the same command strings that the CLI accepts.

### Python:

```python
from gammaloop import GammaLoopAPI

gl = GammaLoopAPI(
    state_folder="./examples/cli/gg_hhh/1L/state",
    boot_commands_path="./examples/cli/gg_hhh/1L/gg_hhh_1L.toml",
)

gl.run("display settings global")
gl.run("run integrate_physical")
```

The Python bindings can be built with:

```bash
just build-api
```

The Python package source is under `crates/gammaloop-api/python/gammaloop`.

There is also a standalone Python API example that builds a differential LU state
and calls `evaluate_sample(...)` / `evaluate_samples(...)`:

```bash
source .venv/bin/activate
NO_SYMBOLICA_OEM_LICENSE=1 SYMBOLICA_LICENSE=... just build-api
cd examples/api/python/epem_a_ddxg_xs_LO
SYMBOLICA_LICENSE=... python inspect_events.py
```

That example uses:
- `examples/api/python/epem_a_ddxg_xs_LO/run.toml`
- `examples/api/python/epem_a_ddxg_xs_LO/inspect_events.py`

and creates its state under:

```text
examples/api/python/epem_a_ddxg_xs_LO/state
```

### Rust:

```rust
use std::path::PathBuf;

use gammaloop_api::{state::CommandHistory, StateLoadOption};

fn main() -> color_eyre::Result<()> {
    let mut loaded = StateLoadOption {
        state_folder: Some(PathBuf::from("./examples/cli/gg_hhh/1L/state")),
        boot_commands_path: Some(PathBuf::from("./examples/cli/gg_hhh/1L/gg_hhh_1L.toml")),
        ..StateLoadOption::default()
    }
    .load()?;

    let command = CommandHistory::from_raw_string("display settings global")?;
    let mut session = loaded.cli_session();
    let _ = session.execute_command(command)?;

    Ok(())
}
```

In Rust, the API entry point is `StateLoadOption::load()`. After loading, you can either start a CLI-style session with `loaded.cli_session()` or use dedicated command/data structures directly on the loaded state.

There is also a standalone Rust API example, designed to run with
[`rust-script`](https://crates.io/crates/rust-script), that builds the same
differential LU state and prints the formatted rich result returned by
`evaluate_sample(...)` / `evaluate_samples(...)`:

```bash
NO_SYMBOLICA_OEM_LICENSE=1 EXTRA_MACOS_LIBS_FOR_GNU_GCC=T SYMBOLICA_LICENSE=... \
rust-script --debug examples/api/rust/epem_a_ddxg_xs_LO/inspect_events.rs
```

That example uses:
- `examples/api/rust/epem_a_ddxg_xs_LO/run.toml`
- `examples/api/rust/epem_a_ddxg_xs_LO/inspect_events.rs`

and creates its state under:

```text
examples/api/rust/epem_a_ddxg_xs_LO/state
```

Current note on shape and precision:

- `evaluate_sample(...)` returns one sample result plus one observable bundle for
  that single-sample batch.
- `evaluate_samples(...)` returns per-sample evaluation results together with one
  batch-global observable bundle merged over all samples in the batch.
- the default Rust/Python `evaluate_sample(...)` / `evaluate_samples(...)` API
  downcasts returned numeric values to `f64`, even when `use_arb_prec=true`.
- for Rust-only workflows that need to keep the retained precision, use the
  parallel `evaluate_sample_precise(...)` / `evaluate_samples_precise(...)`
  helpers from `gammaloop_api::commands::evaluate_samples`; they return
  precision-tagged results instead of widening the Python API away from its
  `f64` contract.

## Examples and tests

Most historical files under `examples/` are being refreshed to match the current CLI.
The main up-to-date CLI example is:

- `examples/cli/gg_hhh/1L/gg_hhh_1L.toml`

Integration tests and additional fixtures live under:
- `tests/tests/`
- `tests/resources/run_cards/`

## Getting help

Use the built-in help for the current command surface:

```bash
./gammaloop --help
./gammaloop help
./gammaloop help generate
```

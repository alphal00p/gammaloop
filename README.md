<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-dark.svg">
  <img src="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-light.svg" width="300">
</picture>
</div>

*Computation of differential cross-sections using Local Unitarity.*

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

Additional requirements for rendering diagrams:

```bash
cargo install clinnet
cargo install typst
```

If local dependency management is inconvenient, you can instead use the provided Nix flake:

```bash
nix develop
```

and then run the usual build commands from within that development shell.

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
- `run.toml`
- `processes/`

## Python API

The Rust CLI is the primary interface under active iteration.
The Python bindings are still available and can be built with:

```bash
just build-api
```

The Python package source is under `crates/gammaloop-api/python/gammaloop`.

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

#import "../../shared.typ": callout, boundary, product-link, source-link

#let overview = [
= Overview

GammaLoop computes differential collider cross-sections with the Local Unitarity method.
It is an application and research codebase rather than a thin numerical library: process
generation, model import, integrand construction, integration, observables, persistence,
and diagnostics meet in one stateful workflow.

#callout("Primary interface", [
  The command-line interface is the primary user interface. A run card creates a state,
  executes commands, and leaves a reusable working directory. The Rust and Python APIs
  load the same persisted state and expose selected structured operations; they do not define
  a second, independent execution model.
])

== A stateful run

A typical source checkout is built and exercised from the repository root:

```sh
just build-cli-release
./gammaloop ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml
```

The run card imports a model, generates the requested process, executes its command blocks,
and records the resulting state. Resume that directory explicitly for subsequent work:

```sh
./gammaloop -s ./examples/cli/gg_hhh/1L/state \
  run integrate_physical -c "quit -o"
```

The persisted `run.toml` records how to replay the run. The settings files describe global and
default runtime configuration, while `processes/` contains process-specific state. Ordinary
runs create `gammaloop_state/` unless a different state path is supplied.

== Installation and external tools

GammaLoop is currently developed primarily from source. The supported development path is
the repository's Nix shell, or a local Rust toolchain together with `just`, a recent GNU
toolchain, and Python 3.11 or newer when building bindings. FORM 4.2.1 or newer is needed
for analytical integration of integrated UV counterterms. UFO model import also needs the
Python `ufo-model-loader` package.

```sh
nix develop
just build-cli
just build-api
```

Diagram rendering is a separate concern and uses Clinnet and Typst. Building the CLI does
not imply that these drawing tools are installed.

== Related crates

#boundary("How the crates fit together", [
  #product-link("linnet", label: "Linnet") provides the half-edge graph model and graph
  algorithms. #product-link("spenso", label: "Spenso") provides typed tensors, tensor
  structures, and network execution. #product-link("idenso", label: "Idenso") handles symbolic
  tensor identities and algebraic simplification. #product-link("vakint", label: "Vakint")
  matches and evaluates vacuum-integral topologies. GammaLoop combines these components into a
  collider-calculation workflow.
])

For a deeper view of how the command layer, state, integrands, and external components fit
together, see
#source-link("docs/architecture/architecture-current.md", label: "the architecture overview").

== Where to begin

- Use built-in `--help` for the CLI options supported by your installed version.
- Start with the included `gg_hhh` run card for a complete state lifecycle.
- Use the Rust or Python API only when a program needs structured access to loaded state or
  sample-evaluation results.
- Check prerequisites before starting an expensive integration; some workflows require a
  Symbolica license or external tools such as FORM.
]

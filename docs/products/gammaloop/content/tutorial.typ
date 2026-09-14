#import "../../shared.typ": boundary, callout, source-link

#let tutorial = [
= Tutorial

This tutorial builds the GammaLoop command-line interface and uses the included one-loop
`g g -> h h h` run card to create a state that can be inspected and resumed. Compilation and
process generation can take longer than a `--help` smoke test.

For the shortest installation and first-run path, begin with
#link("quickstart/cli/")[Using GammaLoop from the command line]. Then return here to create,
inspect, and resume the larger bundled example.

#callout("What this example establishes", [
  The card selects one fixed-helicity, one-loop amplitude contribution and supplies an
  explicit `1/8` color projector. It demonstrates software state creation, inspection, and
  resumption; it is not an unpolarized, automatically color-averaged cross section. Read the
  #link("physics/")[method and capability boundary] and
  #link("guides/conventions/")[normalization conventions] before interpreting either
  integration block as a physical result.
])

== Prerequisites

Work from a GammaLoop source checkout. The supported development environment is the repository
Nix shell; without Nix you need Rust 1.85 or newer, `just`, a recent GNU toolchain, and the
dependencies described in the root README. The bundled Standard Model used below does not need
the UFO loader. A Symbolica license may be required by your build and use case.

// docs-example: syntax
```sh
nix develop
just build-cli
./gammaloop --help
```

The last command should print the one-shot CLI options and command tree. The `./gammaloop`
wrapper selects the binary built under `target/dev-optim/`.

#callout("Verification scope and cost", [
  The docs harness checks these shell blocks for syntax; it does not build or execute
  GammaLoop. Verify a clean checkout manually in `nix develop`. A cold CLI build and process
  generation are the high-resource path and can take tens of minutes, while reopening a valid
  state is normally much cheaper. The invariant for this tutorial is a resumable state with the
  files listed below, not a completed Monte Carlo integration.
])

== Generate an example process

The run card declares a state folder, settings, and named command blocks. Run its `generate`
block and persist the result on exit:

// docs-example: syntax gammaloop-first-state
```sh
./gammaloop --clean-state \
  ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml \
  run generate -c "quit -o"
```

#callout([`--clean-state` removes the resolved state first], [
  Use it only for this first, reproducible run. Omit it when resuming work you want to keep.
  Its generated #link("reference/cli/commands/gammaloop/#argument-gammaloop-clean-state-349074e9fab4dbbb")[argument reference] records the exact
  valueless-flag semantics.
  The example card resolves its state to `examples/cli/gg_hhh/1L/state` and requests ten
  worker cores; make a private copy of the card before changing either value.
])

The block imports `sm-default.json`, generates the selected one-loop pentagon amplitude
contribution, applies the card's explicit helicity and color choices, builds its integrand,
saves DOT data, and writes the state. A successful run exits without an error and leaves at
least `run.toml`, `global_settings.toml`,
`default_runtime_settings.toml`, and `processes/` in the state folder. `run.toml` records the
commands and settings needed to replay the run; it is not merely a log file.

== Resume the saved state

Load the saved state and ask the active session to display its processes:

// docs-example: syntax
```sh
./gammaloop -s ./examples/cli/gg_hhh/1L/state \
  run -c "display processes; quit -o"
```

Success means the display includes the generated `gg_hhh` process and the command exits while
keeping the same state. You can now run the card's `integrate_euclidean` or
`integrate_physical` block, but those are deliberately beyond the first-success path because
they request a substantial Monte Carlo integration and require a separate normalization and
validation plan.

== Inspect before expensive work

Reopen the state read-only and inspect both the generated structure and the effective integrator
settings before launching either integration block:

// docs-example: syntax
```sh
./gammaloop --read-only-state \
  -s ./examples/cli/gg_hhh/1L/state \
  run -c "display integrand -p gg_hhh -i 1L; display settings process -p gg_hhh -i 1L integrator; quit"
```

`display integrand` reports the generation backend, compiled record size, graph and graph-group
counts, and detailed orientation, loop-momentum-basis, cut, and threshold tables. The settings
command resolves the values attached to this generated integrand, which may differ from the
current defaults. Check these records together: an unexpected graph count points back to process
generation, while an unexpected sample budget, target accuracy, or seed belongs to runtime
settings.

#callout("Read-only inspection is not persistence", [
  `--read-only-state` prevents commands from writing inside the active state directory and
  disables file logging there. It does not make the in-memory session immutable. This command
  intentionally exits without `-o`; change settings in a separate writable session or in the
  run card that reconstructs the state.
])

The exact fields and selectors are kept current in the generated
#link("reference/cli/commands/gammaloop/display/integrand/#command-gammaloop-display-integrand-7515b87f5a3f8b18")[integrand display] and
#link("reference/cli/commands/gammaloop/display/settings/process/#command-gammaloop-display-settings-process-66a6624acc30949a")[process-settings] references.

== Troubleshooting and next steps

- If the wrapper cannot find a binary, rerun `just build-cli` and confirm that
  `target/dev-optim/gammaloop` exists.
- If compilation or process generation exhausts local resources, copy the card and reduce the
  values under `cli_settings.global.n_cores`; do not edit a state midway through generation.
- If startup reports a state fingerprint or settings mismatch, replay the card with a new
  state folder instead of transplanting only its `processes/` directory.
- Use `./gammaloop help generate` to inspect the flags supported by your installed version,
  then continue with the #link("guides/process-generation/")[process-generation guide] and the
  #link("reference/cli/commands/gammaloop/generate/#command-gammaloop-generate-9dcc9f488fe75777")[generate-command reference].

#boundary("Example run card", [
  See the complete
  #source-link("examples/cli/gg_hhh/1L/gg_hhh_1L.toml", label: "gg_hhh run card") before
  adapting the process, state location, or resource settings for your own calculation.
])
]

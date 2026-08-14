#import "../../shared.typ": boundary, callout, source-link

#let tutorial = [
= Tutorial

This tutorial builds the GammaLoop command-line interface and uses the maintained one-loop
`g g -> h h h` run card to create a state that can be inspected and resumed. It is a real
process-generation workflow, so expect compilation to take longer than a `--help` smoke test.

== Prerequisites

Work from a GammaLoop source checkout. The supported development environment is the repository
Nix shell; without Nix you need Rust 1.85 or newer, `just`, a recent GNU toolchain, and the
dependencies described in the root README. The bundled Standard Model used below does not need
the UFO loader. A Symbolica license may be required by your build and use case.

```sh
nix develop
just build-cli
./gammaloop --help
```

The last command should print the one-shot CLI options and command tree. The `./gammaloop`
wrapper selects the binary built under `target/dev-optim/`.

== Generate one maintained process

The run card declares a state folder, settings, and named command blocks. Run its `generate`
block and persist the result on exit:

```sh
./gammaloop --clean-state \
  ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml \
  run generate -c "quit -o"
```

#callout("`--clean-state` removes the resolved state first", [
  Use it only for this first, reproducible run. Omit it when resuming work you want to keep.
  The maintained card resolves its state to `examples/cli/gg_hhh/1L/state` and requests ten
  worker cores; make a private copy of the card before changing either value.
])

The block imports `sm-default.json`, generates the selected one-loop pentagon contribution,
builds its integrand, saves DOT data, and writes the state. A successful run exits without an
error and leaves at least `run.toml`, `global_settings.toml`,
`default_runtime_settings.toml`, and `processes/` in the state folder. `run.toml` is the audit
and replay description; it is not merely a log file.

== Prove that the state is reusable

Load the saved state and ask the active session to display its processes:

```sh
./gammaloop -s ./examples/cli/gg_hhh/1L/state \
  run -c "display processes; quit -o"
```

Success means the display includes the generated `gg_hhh` process and the command exits while
keeping the same state. You can now run the card's `integrate_euclidean` or
`integrate_physical` block, but those are deliberately beyond the first-success path because
they request a substantial Monte Carlo integration.

== Troubleshooting and next steps

- If the wrapper cannot find a binary, rerun `just build-cli` and confirm that
  `target/dev-optim/gammaloop` exists.
- If compilation or process generation exhausts local resources, copy the card and reduce the
  values under `cli_settings.global.n_cores`; do not edit a state midway through generation.
- If startup reports a state fingerprint or settings mismatch, replay the card with a new
  state folder instead of transplanting only its `processes/` directory.
- Use `./gammaloop help generate` for flags from the exact build you are running, then continue
  with the process-generation manual and the generated CLI reference.

#boundary("Ground truth for this walkthrough", [
  The complete card is #source-link("examples/cli/gg_hhh/1L/gg_hhh_1L.toml", label: "the maintained gg_hhh example").
  The executable's parser and generated CLI catalog own option spelling; this tutorial owns the
  sequence and the meaning of the resulting state.
])
]

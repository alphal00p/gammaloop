= BNL Contraction Order Performance Notes

#quote(block: true)[
#strong[Status:] Archived, partially reproducible performance evidence

These measurements were captured on 2026-05-27 and 2026-05-28 at evidence revision
`4fdbf430b29edd24b9e1292c07ab24a151427dc7`. The two input atoms and the
`bnl_evaluator_atom_mwe` example remain, but the original note did not identify the machine,
operating system, Nix closure, or Symbolica revision. Re-run the commands to establish a new
baseline; do not compare current timings directly with these values.
]

Input:

```sh
examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_exact.sym first:1 --alias-scalars 4096
```

The runs below used an isolated target directory and low cargo parallelism to avoid
competing with background cargo checks and to reduce OOM risk:

```sh
env CARGO_TARGET_DIR=/private/tmp/gammaloop-codex-target CARGO_BUILD_JOBS=2 \
  cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- \
  examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_exact.sym first:1 \
  --alias-scalars 4096 --contraction-order ORDER
```

The first `dev-optim` build took `18m 34s`. The table timings below are from the
example's own stage timers, so compile time is excluded.

== Direct Concrete Network

#table(
  columns: 3,
  table.header([*Order*], [*Network parse*], [*Network execute*]),
  [sparse-atom-aware], [65.897583ms], [1.52775ms],
  [atom-aware], [63.773667ms], [1.630291ms],
  [entry-aware], [62.855959ms], [2.201458ms],
  [result-rank-only], [60.812834ms], [6.8265ms],
)

#table(
  columns: 4,
  table.header([*Order*], [*Entries*], [*Terms*], [*Total entry bytes*]),
  [sparse-atom-aware], [569], [928], [374721],
  [atom-aware], [569], [928], [374721],
  [entry-aware], [569], [928], [374721],
  [result-rank-only], [557], [920], [9262914],
)

#table(
  columns: 3,
  table.header([*Order*], [*Largest entry bytes*], [*Result bytes*]),
  [sparse-atom-aware], [157199], [157216],
  [atom-aware], [157199], [157216],
  [entry-aware], [157199], [157216],
  [result-rank-only], [2581208], [2581225],
)

`sparse-atom-aware` and `atom-aware` choose the non-blowing contraction shape for
this term. The old rank-only score still produces a single huge scalar entry:
about `24.7x` more post-execute tensor scalar bytes and `16.4x` larger result
bytes than `sparse-atom-aware`.

== Symbolic Then Concrete Baseline

Command:

```sh
env CARGO_TARGET_DIR=/private/tmp/gammaloop-codex-target CARGO_BUILD_JOBS=2 \
  cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- \
  examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_exact.sym first:1 \
  --alias-scalars 4096 --contraction-order sparse-atom-aware --symbolic-then-concrete
```

#table(
  columns: 5,
  table.header([*stage*], [*elapsed*], [*terms*], [*bytes*], [*details*]),
  [symbolic network parse], [4.722125ms], [-], [-], [nodes=49 edges=91 tensors=22],
  [symbolic network execute], [0.475375ms], [-], [-], [aliases retained],
  [symbolic aliased root], [-], [1], [1754], [aliases=1],
  [concrete network parse], [56.123666ms], [-], [-], [nodes=49 edges=91 tensors=15],
  [concrete network execute], [1.21825ms], [-], [-], [sparse-atom-aware],
  [concrete tensor entries after execute], [-], [928], [374721], [entries=569 max_entry_bytes=157199],
  [concrete aliased root], [-], [1], [157216], [aliases=1],
)

At capture time, the direct sparse-aware order reached the same post-execute scalar-growth
class as the symbolic-then-concrete path for `first:1`: `569` entries, `928`
terms, `374721` total tensor-entry bytes, and `157216` result bytes.

== Notes

- I did not run `first:4` or `all` in this pass. The purpose here was to validate
  the scoring change without deliberately triggering the known four-term memory
  pressure.
- At evidence revision `4fdbf430b29edd24b9e1292c07ab24a151427dc7`, the sparse-aware score was
  the default `MinResultRank` preset. The
  `EvaluatorSettings::tensor_network_contraction_order` setting can select the
  other presets for comparison.

== Direct All Follow-Up

Date: 2026-05-28

Command:

```sh
env CARGO_TARGET_DIR=/private/tmp/gammaloop-codex-target CARGO_BUILD_JOBS=2 RAYON_NUM_THREADS=1 \
  cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- \
  examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_exact.sym all \
  --alias-scalars 4096 --contraction-order sparse-atom-aware
```

`first:4` is not a meaningful selector for this saved exact atom: Symbolica parses
it as one top-level term. Running `all` directly, without expanding first, completed
without memory pressure.

#table(
  columns: 2,
  table.header([*Metric*], [*Value*]),
  [Selection], [`all`],
  [Network parse], [169.562583ms],
  [Network execute], [2.275083ms],
  [After-execute entries], [569],
  [After-execute terms], [928],
  [After-execute bytes], [374721],
  [Largest entry bytes], [157199],
  [Result bytes], [157216],
)

== Unfiltered Pre-Network Follow-Up

Date: 2026-05-28

The exact atom above is already a single factored top-level term, so `all` and
`first:1` select the same expression. The larger four-alias case is represented
by `bnl_integrated_evaluator_atom_unfiltered_pre_network.sym`.

Command:

```sh
env CARGO_TARGET_DIR=/private/tmp/gammaloop-codex-target CARGO_BUILD_JOBS=2 RAYON_NUM_THREADS=1 \
  cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- \
  examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym all \
  --alias-scalars 4096 --contraction-order sparse-atom-aware
```

#table(
  columns: 2,
  table.header([*Metric*], [*Value*]),
  [Input], [unfiltered pre-network],
  [Aliases], [4],
  [Network parse], [173.041791ms],
  [Network execute], [4.331541ms],
  [After-execute entries], [1669],
  [After-execute terms], [2460],
  [After-execute bytes], [1493281],
  [Largest entry bytes], [545943],
  [Result bytes], [545960],
)

The unfiltered case completes without resolving aliases and without the previous
runaway scalar growth: the largest final entry is about `546kB`, and the aliased
root stays at `546kB` with four retained scalar aliases.

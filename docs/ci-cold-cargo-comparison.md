# Cold main runtime tests: Cargo and Nix

On 11 September 2026, the same main application code completed the selected CI
runtime tests in **21m46s with direct Cargo** and **28m48s through Nix**. Nix took
**7m01s longer (32.3%)** in this single local pair. Both built their Python extension
inside the timed workload and passed the same **1,633 tests**, with **155 skipped**.

| Cold compilation and runtime tests | Direct Cargo | Current Nix implementation |
|---|---:|---:|
| Completion | 1,306.21s / **21m46s** | 1,727.60s / **28m48s** |
| Tests passed | 1,633 | 1,633 |
| Skipped cases | 155 | 155 |
| Cargo compilation during test execution | 0 | 0 |

This covers the seven main runtime groups: core, integration, Python API, Clinnet,
Linnet, Spenso and Vakint. It excludes Clippy, doctests, formatting, graph checks,
release packaging and WASM. The earlier 38m21s Nix measurement included additional
checks and is a separate experiment.

## What was measured

Both routes use the application and CI configuration at
`e56b04d650895fe78f564fea9c4e8dd9b0f3dfa2`, on the main application layout. Initial
source-file hashes match exactly. The Nix runtime roots are identical to those in
the preceding validated candidate; only the selected benchmark scope changes.

- Each route receives the same eight CPUs, 40–47. Cargo uses eight build jobs;
  Nix allows four builders, each configured for eight cores, within that same CPU
  allocation. A running Nix compiler's affinity was verified directly.
- Rust 1.97.0, Cargo 1.97.0, nextest 0.9.140, and the `ci-optim`/`ci_gammaloop`
  profiles are identical. Cargo enables incremental compilation, and Nix package
  producers also generate compiler state. Neither starts with any to reuse.
- Cargo starts with an absent target directory. Nix starts in a new disposable
  store, with every project output and compiler-state output verified absent.
  Both compile Rust dependencies as well as workspace crates.
- Toolchains, native libraries and downloaded crate sources are already available.
  Cargo runs offline; Nix has no substituters, remote builders or upload hook.
- Direct Cargo uses the pinned development environment's compiler and libraries,
  with a normal mutable target directory. It has no Nix build or archive steps.

The runner is Cargo/nextest, preserving the project's CI filtering and scheduling
rules. A literal `cargo test` invocation has a different test runner and, without
`--workspace`, selects only the default member, `gammaloop-api`.

The native Python build uses the Nix module's exact package/feature arguments,
stable Python ABI and interpreter. The subsequent workspace test build uses the
CI test interpreter and enables `python-api-tests` and Symbolica's CI tracing
feature. All required CI features are present. Workspace feature unification and
Nix's package contexts still produce different Cargo compilation-unit graphs.
All 1,788 listed cases match by package, binary, target kind, name, test kind,
ignored flag and selection; 1,633 are selected in each route.

## Where the time went

Cargo runs these stages sequentially, sharing one initially empty target:

| Cargo stage | Time |
|---|---:|
| Build the Python extension and its dependencies | 7m31s |
| Package the extension | 0.13s |
| Compile the workspace tests, reusing the Python build where possible | 11m15s |
| Execute all selected tests, including Python | 3m00s |

An additional 1.13s inventory-validation step is excluded from the comparison.
Including it gives the initially reported 21m47s. No preliminary `cargo check`
warmed the timed target; the unchanged source had already passed the earlier
formatting, Cargo check and Clippy validation.

Nix overlaps compilation and runtime groups. Clinnet finishes at 10m08s, Vakint at
11m26s, Linnet at 14m41s, Spenso at 16m56s, core at 24m23s and Python at 27m45s.
Integration starts at 25m29s, runs for 3m18s and finishes last, at 28m48s.
These overlapping durations must not be added to obtain completion time.

Nix's final Python packaging takes 8.79s with zero Cargo compilation. Its test
archive steps and runtime checks also compile nothing. The extra cold time occurs
before and alongside test execution; it cannot be attributed to NixCI transfers
or queues because neither is present in this experiment.

Nix evaluation takes 12.89s and copying external inputs/derivation definitions
takes 9.12s, separately from the timed build. An exact Nix repeat takes 0.23s and
executes no builders or tests; it reuses successful results.

## Evidence and limits

The compact [JSON results](ci-cold-cargo-comparison.json) retain timings, coverage,
initial-state checks and commands. Raw plans, frozen sources, logs, inventories,
the disposable store and `comparison.csv` are in `/tmp/ci-cargo-nix-cold-tests`.

This is one sequential pair on a shared host. CPU allocation is matched, but
memory/storage contention and other host activity are not isolated; operating
system file caches were not flushed. Nix's summed builder occupancy is not CPU
time or a cost estimate. These measurements establish a matched local cold
baseline, not a universal speed ratio. No NixCI suite was pushed for this test.

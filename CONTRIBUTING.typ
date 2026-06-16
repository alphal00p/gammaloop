= Contributing

These guidelines are shared by human contributors and coding agents. They cover
the project workflow, the engineering style we want, and the subsystem-specific
constraints that are easy to miss during refactors.

== Core Engineering Principles

These rules are intentionally broad and should shape most code changes.

- Prefer one semantic boundary per concept. Before adding a trait, helper, mode,
  or adapter, look for the existing abstraction that already owns the behavior
  and extend that instead. Before adding helpers, structs, or methods, search
  similar use cases; ask maintainers only if ownership or duplication remains
  unclear after that search.
- Collapse duplicate abstractions aggressively. If two names describe the same
  responsibility, merge them and update call sites rather than keeping
  compatibility shims.
- Prefer methods, trait methods, or impl-associated functions over bare free
  functions. A helper should live on the type or trait whose invariant it uses.
- Do not preserve internal API compatibility unless explicitly requested.
  Breaking internal Rust/Python APIs is acceptable when it removes duplication
  or clarifies ownership.
- Avoid compatibility fallbacks and "old path plus new path" designs. Pick the
  clearer model and migrate usages.
- Keep changes minimal but complete: remove obsolete code paths, imports, docs,
  and call sites in the same change. Preserve useful comments and rationale
  when moving code; remove demonstrably obsolete comments with the behavior
  they describe, and ask if their relevance is unclear.
- Prefer concise, idiomatic code and remove duplication. Review the diff for
  readability before finishing; line counts are not a quality target.
- Prefer concrete return types over sentinel-style APIs. Avoid `Option` where
  `None` means "nothing happened" if the caller always needs a usable result.
- Make transformations explicit and idempotent. If recursive rewriting is
  involved, include a fixed-point check or another clear loop guard.
- Avoid relying on fragile errors for control flow when a direct predicate or
  structural check would express the intent better.
- Do not add broad helper layers just to avoid touching call sites. Update the
  call sites when the API should change.


== Investigation Discipline

- Before fixing a suspected bug, find the existing code path that should own the
  behavior.
- When asked to identify failure modes, do not patch yet; report the concrete
  failing cases and where they enter the code.
- When asked for a plan, do not make code changes until the plan is accepted.
- If a requested change looks like only a rename, stop and check whether there
  is a deeper design change implied.
- Ignore unrelated untracked local artifacts by default (editor swap files,
  scratch docs, local example edits, profiling outputs, etc.) unless the task
  clearly requires them.

=== Discrepancy Triage
<discrepancy-triage>
Before investigating an implementation mismatch at a low level, first
perform a bird\'s-eye comparison of the two complete pipelines. This
step is mandatory: it prevents a difference in a shallow representation
boundary from being debugged as a problem in shared algebra or physics.

+ Write down the inputs, requested modes, and ordered processing stages
  for both routes. Mark every stage as shared or different, and identify
  the last shared boundary and the first boundary that can produce
  different state.
+ Record the smallest result matrix that distinguishes the routes.
  Toggle one independent feature at a time, such as empty versus
  nonempty UV forest, local versus integrated counterterms, threshold
  subtraction, cut, residue order, LMB channel, orientation, evaluator
  mode, and numeric precision.
+ Apply logical exclusions before proposing causes. A mechanism shared
  by two routes cannot explain a difference between them unless the
  inputs reaching that mechanism have already diverged. An independently
  validated downstream engine is not a candidate until its two actual
  inputs are shown to differ or an identical-input A/B test fails.
+ Inspect the representation at the first differing boundary before
  tracing deeper code. Prefer durable artifacts such as
  `save standalone --json`, UV forest exports, generation reports, and
  structured debug logs. Compare branch counts, keys, maps, selectors,
  and factorized coefficients before comparing only final floating-point
  totals.
+ Reduce the mismatch hierarchically: total, graph, forest, cut, residue
  order, LMB, orientation, then individual term. Stop as soon as the
  first unequal pair is found and make that pair the reproducer.
+ For selector-local versus explicit-sum representations, compare
  `sum(selector * body)` with `sum(body)` directly. Evaluate the
  complete selector truth table and prove that each explicit branch is
  selected exactly once. Do not investigate contour, residue,
  reconstruction, or numerator machinery until this shallow
  partition-of-unity comparison passes.

Diagnostic normalization must preserve the factorization of graph
numerators, including test-only copies. Denominator-only algebra and
finite tensor-component contraction are distinct from distributing a
graph numerator. Each progress report for a discrepancy should state:
what is shared, what first differs, what has been excluded, the smallest
current reproducer, and the single next comparison that will reduce it.

==== Projected local-4D UV reconstruction stop rule
<projected-local-4d-uv-reconstruction-stop-rule>
For a projected local-4D UV-to-CFF mismatch, reconstruction must be
certified #emph[before] investigating CFF recursion, contour signs, or
residue aggregation. This is a mandatory correctness boundary, not an
optional diagnostic:

+ Retain raw post-Taylor sectors in the compatible hard sub-LMB,
  including original edge owners and provenance roles for subsequent
  outer Taylor operations. Build canonical algebra only in a separate
  projection view.
+ Resolve denominator classes by exact signed routing, mass, full
  polynomial (including its prescription), and component domain. Keep
  physical source occurrences as a separate incidence witness. Construct
  the UV skeleton from that witness; derivatives may add serial copies
  of retained lines. Contract equivalent channels only with a certified
  serial-path or pure-cycle incidence. Nonadjacent equivalent channels
  retain their original incidence; do not infer a graph from a signature
  matrix.
+ Completed hard provenance roles zero and one may transfer to a
  certified denominator class and use occurrences of that class,
  including original numerator factors. Preserve the literal
  odd-momentum sign independently of denominator evenness. A hard
  carrier without an unambiguous surviving pole class stays a fixed
  affine carrier; physical-source and soft provenance keep their
  existing meanings. New soft factors retain explicit provenance until
  outer-CFF assembly; only those factors may use an exactly certified
  off-shell routing through active cograph edges, including every fixed
  external shift.
+ Apply the exact immutable production assignment plan, including the
  signed hard/raw/parsed conversion `H = h R`, `P = r R`, and
  `H^0 = h r P^0`.
+ Substitute the source post-Taylor numerator and the reconstructed
  UV-EMR numerator into one neutral set of formal loop four-momenta (and
  the same fixed external data) and require their exact symbolic
  difference to vanish.
+ Independently require equality of denominator momentum, mass,
  multiplicity, and component domain, using `D(Q) = D(-Q)` only on the
  denominator side.

Only after both exact certificates pass may a discrepancy be attributed
to generalized-CFF input normalization, CFF generation, component
composition, or residue aggregation. The common LMB is only a coordinate
chart for this proof; it never supplies EMR ownership or CFF rank
capacity. Correctness of the EMR rewrite precedes dispatch optimization.
Production compares at most three certified proposals by their actual
native generated source-map row count, before surface conversion or
cut/host selection. Rank orders the proposals and breaks count ties; it
is not a substitute for that count or a proof of global optimality. The
selected payload and its exact assignment must stay together. Graph
numerators remain factorized, including diagnostic copies. Establish the
identity with factor-preserving rewrites and exact cancellation;
evaluations at explicit coordinates provide additional checks, not a
general symbolic proof. The worked GL04 `1zs/T2` certificate is
maintained in
#link("docs/architecture/exact-powered-denominator-cff-lifting.typ#worked-live-reproducer-gl04-temporal-square-1zst2")[`docs/architecture/exact-powered-denominator-cff-lifting.typ`];.

== Debug Logging Pattern

- Prefer `debug_tags!` plus the log filter environment variables over ad hoc
  debug files, one-off `println!`, or custom dump env vars. This keeps
  investigative output routable to display and JSONL log sinks with the same
  filtering syntax.
- Add semantic tags that describe the pipeline and subsystem, for example:
  ```rust
  debug_tags!(#uv, #integrated, #vakint, #trace;
      stage = "to_vakint_integrand_after_den_to_prop",
      reduced = %reduced_label,
      display.preview = %expr.log_print(Some(120)),
      file.expr = %expr.to_plain_string(),
      "Vakint trace"
  );
  ```
- Use a `stage` field for fine-grained ordering within a code path. Use
  concise display fields for terminal output and `file.*` fields for large or
  exact payloads. `file.*` fields are hidden from the display sink and written
  as structured JSONL logfile fields with the routing prefix stripped
  (`file.expr` is stored as `expr`).
- To capture a focused test trace to file while keeping stderr quiet, run for
  example:
  ```bash
  GL_DISPLAY_FILTER=off \
  GL_LOGFILE_FILTER='[{#uv,#integrated,#vakint,#trace}]=debug' \
  GL_TEST_LOG_DIR=/tmp/gammaloop-test-logs \
  cargo test -p gammalooprs --test test_renormalization --profile dev-optim \
    finite_part_ghost_3loop_renormalization_sum_origin_mwe -- --ignored --nocapture
  ```
- Use `GL_ALL_LOG_FILTER` only when display and logfile should receive the same
  filter. Otherwise prefer `GL_DISPLAY_FILTER` and `GL_LOGFILE_FILTER` so large
  `file.*` payloads can be kept out of terminal output.

== Code Quality Baseline

=== Rust

- Format Rust changes and run checks appropriate to the affected code. Use
  `cargo check` for quick feedback when useful; it is not a prerequisite for
  every build or test run.
- Run Clippy for relevant Rust changes and address warnings where practical.
- Naming: `snake_case` for functions/modules, `CamelCase` for types/traits.
- Prefer narrow imports and method-call syntax over repeated fully qualified
  paths; retain qualification for disambiguation or macro hygiene.
- Prefer methods on types to bare functions.
- Prefer explicit function composition over function nesting: when a value
  crosses a semantic boundary, bind it or use `From`/`TryFrom` so the
  conversion is visible at the call site.
- Name methods for the semantic result they compute, not for the current caller.
  Avoid `*_for_<caller>` methods unless the result has a genuinely
  caller-specific invariant.
- For strategy/configuration variants, prefer a small builder-style API on the
  owning type over extra method parameters or `*_with_*` methods. Named
  constructors are enough when there are only a few modes; use a fuller builder
  only when options become independent or numerous.
- Keep shared behavior behind the primary trait method when a type implements a
  trait. Put mode selection in the implementing type's configuration rather
  than exposing parallel helper methods that bypass the trait.
- Keep mode and strategy enums private unless callers need to inspect or match
  on them.
- Follow clippy suggestions closely.
- Prefer `pub(crate)` unless fully public exposure is necessary; this keeps
  unused components visible to clippy.

=== Python

- Use the ruff formatter.

=== Tests

- Never weaken tests merely to make failures disappear. Updating expectations
  for an explicitly requested behavior change is permitted; otherwise ask
  before changing what a failing test asserts.
- During development, run checks relevant to the change. Before final review,
  run the selected full CI suite as described below.
- Install `cargo-nextest` 0.9.115 or newer. The repository configuration
  enforces this minimum; update an existing installation with
  `cargo nextest self update`.
- Rust integration tests live in `tests/` with shared fixtures in
  `tests/resources/`.
- Use `cargo nextest` via `just test TEST_NAME` for targeted integration runs.
- Pass emulated libtest arguments after `--`, for example
  `just test_gammaloop -- --skip TEST_NAME`.
- Put targeted unit tests next to the implementation they exercise.
- Use broader integration tests only for cross-module behavior.
- Add tests for edge cases that motivated the change, especially when collapsing
  duplicated logic.
- Numerator oracles must retain physical locality: edge factors depend only on
  their edge momentum, and vertex factors only on incident momenta. Non-local
  algebraic diagnostics do not certify UV subtraction or graph reconstruction.

=== Docs

- Documentation should explain conventions and algorithmic intent, not restate
  function names. Avoid tautological docs like "apply the convention"; say what
  the convention does.
- Keep docs synchronized with API refactors; remove references to deleted
  traits, functions, and modes.

==== Documentation Source Format

- Typst (`.typ`) is the canonical source format for source-controlled authored
  repository documentation, including contributor guidance, manuals,
  architecture notes, investigations, plans, and changelogs. Never add a new
  Markdown/HTML prose source or a parallel copy.
- Files named exactly `AGENTS.md` and `README.md` are the only Markdown
  compatibility exceptions: agent tooling discovers the former, while GitHub,
  Cargo, and package indexes discover the latter. Keep those files concise and
  link them to the canonical Typst documentation.
- Do not add a separately maintained `README.typ`; the compatibility README is
  its own deliberately small source. This applies in every directory. Changing
  extension spelling or case (`.MD`, `.mdx`, `.HTML`, `.TYP`, and similar) does
  not create another exception; canonical Typst sources use lowercase `.typ`.
- Existing Markdown and HTML prose is legacy migration debt, not precedent for
  new files. Any change to its authored content must migrate the canonical
  source to Typst and update every registry and link in the same change instead
  of retaining both formats.
- Build outputs may use a downstream-required format, but generated output is
  not a second source and must not be checked in as hand-maintained prose.
- Architecture docs live in `docs/architecture/` and are classified by
  `docs/developers.toml`. New entries must point directly to `.typ` sources.
  Each source begins with exactly one `= Title`, matching its registry title;
  the builder supplies the page-level heading and verifies that contract.
- `docs/legacy-prose.toml` records that the pre-policy prose migration is
  complete. Its inventory must remain empty; the documentation check rejects
  additions and parallel editable copies.

==== Typst Graph Debugging

When debugging Linnest layout output from Typst, expose the graph records through
Typst metadata and inspect them with `typst query`. The metadata element must be
wrapped in markup, because Typst labels can only be attached in markup mode.

```typst
[#metadata(graph.nodes(g)) <linnest-tree-nodes>]
[#metadata(graph.edges(g)) <linnest-tree-edges>]
```

Then query the values from the same root used for compilation:

```bash
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-nodes>' --field value --pretty
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-edges>' --field value --pretty
```

This is useful for checking generated `pos`, `label-pos`, `layout-width`,
`layout-height`, `label-width`, and `label-height` values without inferring them
from the rendered PDF.

== Repository Map

- `crates/gammalooprs/src/` holds the core Rust implementation.
- `crates/gammaloop-api/` is the workspace member for the Rust API and Python
  bindings (`crates/gammaloop-api/src` and
  `crates/gammaloop-api/python/gammaloop`).
- `tests/` contains the integration-test crate and shared fixtures.
- `examples/` includes command cards and notebooks for end-to-end runs.
- `bin/` provides scripts like `compile.sh`, `run_tests.sh`, and the
  `gammaloop` entrypoints after build.
- `assets/` stores schemas/data and model files (`assets/models/`).
- Rust benches live under `crates/gammalooprs/benches/`.

== Common Commands

- Build Rust CLI in dev-optim:
  ```bash
  just build-cli
  ```
- Build Rust CLI with stable Python ABI:
  ```bash
  just build-cli-abi
  ```
- Build Python API bindings:
  ```bash
  just build-api
  ```
- Build Python API bindings with stable ABI:
  ```bash
  just build-api-abi
  ```
- Format and lint:
  ```bash
  just fmt
  just clippy
  ```
- If you run into a macOS linking issue complaining about a missing `__emul...`
  symbol, try building with `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T`; see `build.rs`
  for the impact of this setting.

The root `justfile` keeps build and lint commands. Test, NixCI and drawing recipes
are imported from `just/tests.just`, `just/ci.just` and `just/drawing.just`; run all
commands from the repository as before. Use `just --list` to see them.

== CI Readiness
<ci-readiness>

The `final-review` label identifies PRs ready for final review. Merging requires
an open, non-draft PR with this label, explicit top-level `enable = true` in
`nix-ci.nix`, and successful required CI checks. Draft and unlabeled PRs are
valid work in progress; their readiness gate intentionally fails to block merging.

During implementation, agents set `enable = false` unless explicitly instructed
otherwise. Human contributors may enable CI earlier. Read-only investigations
do not change the toggle. Keep `enable = true` on `main`.

Before requesting final review:

+ Set top-level `enable = true` in `nix-ci.nix`.
+ On `itphlies`, reuse/download matching cache outputs through Nix substitution.
  Run the selected full CI suite against the final code and upload successful
  results when credentials are available, following the cache workflow below.
+ Push the validated code, mark the PR non-draft, then apply `final-review`.

When returning to development, remove `final-review` and commit
`enable = false`. Removing the label or returning to draft cancels superseded
runs of the Nix and Continuous integration Actions workflows. It does not cancel
already-running NixCI jobs: the committed toggle controls subsequent NixCI work.
Those two Actions workflows run automatically for non-draft, labeled PRs to
`main`; main pushes, merge-group runs, and manual dispatch remain available.
Other workflows retain their own triggers.

Edit only the top-level toggle manually: `just ci-update` preserves it while
regenerating the scheduling configuration. Missing or nonboolean toggles and
other generated-content drift fail validation. `dependency-discovery.enable`
is independent and must not be used to disable CI. Local checks and uploads
remain available with NixCI disabled.

== NixCI Cache
<nixci-cache>

See #link("docs/architecture/ci.typ")[CI maintenance and measurements] for the build design, measured
results, branch migration and `just ci-report` usage.

Before final review, contributors and agents must run the selected CI checks
locally and upload successful results when cache credentials are available.
Ordinary development pushes require checks relevant to the change, not the full
suite and upload cycle.
Work from the repository root; enter `nix develop` if you need the pinned Just
and build tools. Stage newly added files first so the Git-backed flake sees them.

`just check` runs Cargo checking only. For the selected CI suite without uploading:

```sh
just ci-checks
```

This runs the selected Rust and Python tests, Clippy, doctests, formatting and
workspace/CI graph checks. Licensed checks require `SYMBOLICA_LICENSE` in your
local environment. An unchanged successful Nix check can be reused without
executing its tests again. `nix flake check --impure` additionally builds the CLI,
documentation and WASM checks exported by the flake.

On `itphlies`, run the following before pushing CI-enabled work, including
updates to a PR already in final review. The host's Nix daemon is configured to
download from the NixCI cache using its signing key and machine credentials.
Leave substitution enabled: the Nix builds in this command automatically reuse
local outputs or download matching cached outputs, build/check what is missing,
then upload the selected successful outputs and producer closures. There is no
need to download the entire cache or force rebuilds. Completing this locally
before the push minimizes work left for NixCI.

Use the same command on other hosts for the final-review pre-push workflow when
cache credentials are available:

```sh
just ci-checks-and-upload
```

This already runs `just ci-checks`, so there is no need to run both commands.
Wait for the command to succeed and print `CI upload completed` before pushing.
Run it against the final code you intend to push; rerun it if that code changes.
If you lack upload credentials, run `just ci-checks` and mention in your handoff
or PR that the results were not uploaded.

After changing dependencies or features, run `cargo hakari generate` and
`cargo hakari verify` in `nix develop`. Commit the updated workspace-hack manifest
and lockfile. After changing workspace manifests or test groups, run
`just ci-update` before the checks and commit the regenerated
`nix/ci-workspace-graph.json` and `nix-ci.nix` with the change.

Put your NixCI token in `~/.netrc`, or point `NIXCI_NETRC` at an existing file:

```sh
export NIXCI_NETRC=/path/to/existing/netrc
```

On Clan-managed hosts, use the configured `NIXCI_NETRC`; a fresh login may be
needed after upload access is enabled. Keep tokens and license values outside
Git. See #link("https://nix-ci.com/documentation/nix-ci-cache")[the NixCI cache documentation].
Uploads run as your user and need no additional Nix privileges. Entering the dev
shell enables no upload hook. Bare Cargo builds do not populate this cache.
An untrusted daemon may warn that it ignores `netrc-file`; the upload client
still uses that file to authenticate with the destination cache.

The upload command first requires successful checks, then realizes the outputs
selected by `nix/ci.nix` and publishes their runtime closures. These explicit
producer targets retain the binaries and compiler artifacts that test-result
outputs alone would omit. Existing local outputs are published too. It does not
select packaging, documentation, WASM, dev shells or unrelated repositories.
Required shared dependencies can still be part of the uploaded closures.

Check time, publication preparation and upload time are reported separately.
An upload failure returns a nonzero status; fix the reported problem and rerun
`just ci-checks-and-upload`. Successful checks remain cached, so retrying does not
require rerunning matching tests. Use plain `just ci-checks` for local benchmarks
without publishing. Remote reuse also requires matching outputs: a macOS build,
for example, does not replace a Linux build.

The flake retains NixCI's substituter and signing-key settings for downloads.
On a shared daemon, an administrator must configure the cache's trust and access
credentials for substitution; enabling user uploads does not configure downloads.

== CI Completion
<ci-completion>

Prefer completion notifications over an agent repeatedly polling local checks,
cache transfers, or remote CI. When the execution environment supports it, run
long commands with captured logs and an exit-status notification, and resume
only on completion, failure, or required input. Local checks and uploads must
still succeed before the dependent push; background execution is not evidence
of success.

After a push, record the commit and CI run or PR URL, then yield instead of
actively tracking unchanged status unless explicitly asked to monitor. Use an
available event-driven completion wake-up to resume the agent. Verify that the
notification is actually registered; if no such mechanism is available, state
that limitation in the handoff instead of promising an automatic follow-up.
A scheduled monitor still polls; do not substitute it for a completion event
without making that distinction explicit. On resuming, verify that the result
belongs to the intended commit before reporting success or investigating failure.

== Version Control Workflow

=== jj

- If the repo is colocated with jj metadata, prefer `jj` over `git`.
- Start new tasks with `jj new` and immediately
  `jj describe -m "<plan>"`.
- Use `jj log --no-pager` and `jj diff --no-pager` to avoid paging.
- Check the current change before editing; if the working copy has changes but
  no description, run `jj diff` then `jj describe -m "<description>"`.
- Update the current change description with `jj describe -m "<description>"`
  as the plan evolves.
- Use `jj diff` to see changed files.
- Reference: #link("https://jj-vcs.github.io/jj/prerelease/cli-reference/")[Jujutsu CLI reference].

=== Commits And PRs

- Use descriptive branch names without agent prefixes, for example
  `ci-final-review-readiness`. Only use `codex/` when explicitly requested.
- Commit messages are short and descriptive, typically lowercase without scopes
  (for example, `remove edge quotes`).
- For user-requested changes, describe the requested intent and resulting
  behavior; do not automatically copy the full user prompt.
- If there is already a commit, prefer making a new one.
- PRs should include a clear summary, test command(s) run, and any example
  command card used to validate behavior.
- If changes touch CLI state, mention `gammaloop_state/` impacts or migration
  steps.

== Internal API Policy

- API stability is currently a non-goal for internal development: prioritize
  maintainability and structure over preserving external-facing module paths.
- Breaking Rust/Python API changes are acceptable when they simplify
  architecture; document major moves in the change description.
- Backward compatibility with old on-disk GammaLoop states is also currently a
  non-goal: it is acceptable to rename/remove legacy state files, layouts, and
  loaders when that simplifies the implementation.
- Do not keep compatibility fallbacks for old state formats/names unless the
  task explicitly asks for migration support.
- Do not hesitate to modify existing functions when that produces a cleaner and
  more concise design than preserving compatibility or leaving old code
  untouched.

== Project-Specific Notes

The following sections are narrower invariants. They are not general style
preferences; they protect behavior that has been easy to regress.

=== Numeric Precision

- In generic or arbitrary-precision code, do not introduce constants or
  intermediate values through `from_f64`, `std::f64::consts::*`, or other lossy
  `f64` routes unless that exact location is an explicit `f64` boundary by
  design, such as persisted settings that are already `f64` or
  histogram/output accumulation that is intentionally `f64`.
- `f64` values originating from user-supplied settings are an allowed boundary:
  converting those setting values into `F<T>` with `from_f64` is acceptable. Do
  not extend that exception to internally generated constants or intermediate
  values.
- When working with `F<T>` or other precision-generic numeric code, build
  constants from an in-scope representative value of the correct type using
  helpers such as `.zero()`, `.one()`, `.epsilon()`, `.PI()`, `.TAU()`,
  `.from_usize()`, `.from_isize()`, etc., so the active precision is preserved
  exactly.
- If a needed constant or operation cannot be expressed without a dubious
  precision-losing conversion, stop and ask for clarification instead of
  guessing.
- Do not explicitly bring `symbolica::domains::float::FloatLike` into scope in
  GammaLoop code. Prefer the GammaLoop wrappers in `crate::utils`, and if
  `F<T>` is missing an operation needed for precision-safe generic code, add or
  use a helper there instead. If it genuinely looks unavoidable to use
  Symbolica's trait directly, stop and ask first.

=== Settings Serialization

- The settings serialization model relies on the `SHOWDEFAULTS` escape hatch in
  `gammalooprs::utils::serde_utils` when writing persisted state defaults and
  when building completion catalogs from serialized settings.
- When adding or changing settings fields, verify that `Default::default()` and
  the corresponding `skip_serializing_if` helper stay aligned.
- Custom helpers must route through `show_defaults_helper(...)` or
  `IsDefault::is_default`; otherwise the field will disappear from saved
  `global_settings.toml` / `default_runtime_settings.toml` and from completion.

=== State And CLI

- Treat saved-state detection as manifest-based only. Empty folders or folders
  containing only transient runtime artifacts such as `logs/` are scratch state,
  not legacy state, and must be treated as blank state.
- Do not introduce writes into the state folder outside explicit `save state` /
  `quit -o`, except for logfile tracing when that logger is actually enabled.
- Honor `--read-only-state` consistently. When it is enabled, do not write into
  the state folder and prefer cwd-based fallbacks for transient artifacts that
  would otherwise default into the state.
- The run-card field is `command_blocks` with singular `command`, not
  `commands_blocks`. No compatibility alias is kept.
- Integration workspace resume must restore the exact state at the end of the
  last completed iteration. The authoritative resume data lives in the
  workspace manifest/internal state plus per-slot settings and internal
  observable checkpoints, not in user-facing `integration_result.json`.
- Treat workspace snapshots as completed-iteration checkpoints only: do not
  persist partial-iteration state, and keep user-facing latest snapshots plus
  optional `_iter_0001` history separate from the internal resume state.
- CLI runs create a `gammaloop_state/` directory by default; keep it out of
  commits unless intentionally sharing a reproducible state.
- Use `./gammaloop -s <state_dir>` for isolated experiments.

=== Differential LU / Event Processing

- The GammaLoop current-architecture record registered in
  `docs/developers.toml` is the implemented summary for the differential LU
  stack. Keep it synchronized when changing selectors, observables, event
  grouping, or sample-evaluation output.
- Event grouping semantics are by graph-group, not just by graph: if multiple
  LU graphs share the same `group_id`, all accepted cuts from all of those
  graphs belong in the same retained `EventGroup`.
- Histogram error propagation is based on Monte Carlo samples, not on individual
  observable entries. Treat the full retained event-group list from one
  evaluation as one statistical sample for each histogram.
- Histogram accumulation is sparse. Do not explicitly zero-fill untouched bins.
  The sparse state is:
  - histogram-level `sample_count`;
  - per-bin raw stats `entry_count`, `sum_weights`, `sum_weights_squared`,
    `mitigated_fill_count`;
  - underflow and overflow as full bins with the same raw-stat structure.
- Histogram snapshots are intentionally raw-stat snapshots, not presentation-only
  views. They must remain mergeable and reconstructible back into live histogram
  state. Rust-side helpers already exist for `merge(...)`, `merge_in_place(...)`,
  `rebin(...)`, and reconstruction into accumulator state.
- `general.generate_events` controls whether events are retained in returned
  results; observables and selectors may still force temporary event generation
  internally when this is `false`.
- `general.store_additional_weights_in_event` controls whether
  `additional_weights` are stored on events. The map is a `BTreeMap` keyed by
  lightweight identifiers such as `Original`,
  `ThresholdCounterterm { subset_index }`, and `FullMultiplicativeFactor`.

=== API / Python Interop

- `python_api`: Basic Python API support, enabled by default for
  `gammaloop-api`.
- `python_abi`: Enables stable Python ABI (abi3) for cross-version
  compatibility.
  - Use for distribution builds that need to work across Python versions.
  - Adds some performance overhead compared to version-specific builds.
  - Not enabled by default; opt in with `--features python_abi`.
- The Rust API has precise endpoints `evaluate_sample_precise` and
  `evaluate_samples_precise`; the Python API intentionally stays `f64`-only.
- Python API end-to-end tests are subprocess-based and assume the contributor
  already built/installed the Python extension in the active environment
  (`maturin develop` or `just build-api`). Do not try to embed Python into the
  Rust test binary for these tests.
- The top-level `.venv` in the repo is a reasonable default development
  environment for running the Python API examples/tests after building the
  extension there.
- If the extension module environment is an issue when compiling the Python API,
  try `build-api-abi`.

=== Environment And Build Quirks

- The `gammaloop` CLI activates its crate-bound application (OEM) license with
  `symbolica::set_application_key!` before calling the reusable API.
  `NO_SYMBOLICA_OEM_LICENSE` is read there via `option_env!`, so
  opting out requires setting it when compiling the binary, not at runtime.
  Rust/Python libraries and their tests do not activate this OEM license; use a
  user license through `SYMBOLICA_LICENSE` for those entry points.

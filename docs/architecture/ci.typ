= CI maintenance and measurements

Use `just ci-checks` for the selected local suite and `just ci-checks-and-upload`
before final review. The latter runs static checks and publishes compiled runtime
producers without running the runtime suite. NixCI runs that suite after the push.
See #link("../../CONTRIBUTING.typ#nixci-cache")[CONTRIBUTING.typ] for credentials,
licenses and retries. `just check` runs Cargo checking only.

On `itphlies`, run `just ci-checks-and-upload` before pushing CI-enabled work.
The configured Nix daemon downloads matching cached outputs as needed; the
command checks static inputs, builds missing runtime producers and uploads their
closures. Finish the upload before pushing so NixCI can reuse that
work. Keep substitution enabled rather than downloading the whole cache or
forcing rebuilds.

== Configuration and cache boundaries

#table(
  columns: (1fr, 1fr),
  table.header([File], [Responsibility]),
  [`flake.nix`], [Compose public packages, checks, apps and development shells],
  [`nix/rust-workspace.nix`], [Rust sources, features, build profiles and artifact reuse],
  [`nix/documentation.nix`], [Documentation sources, reusable build artifacts and publication checks],
  [`nix/ci.nix`], [Test groups, required producers and NixCI scheduling],
  [`nix/ci-workspace-graph.json`, `nix-ci.nix`], [Generated workspace graph and self-contained NixCI configuration],
  [`just/ci.just`], [Local checks, uploads, regeneration and reporting],
  [`.github/scripts/ci-report.mjs`], [Collect and compare existing CI runs],
)

Automatic NixCI work covers tests, Clippy, doctests, formatting, graph validation,
the full five-product documentation check, necessary producers and final success.
The documentation check validates generated references and examples, renders
latest and snapshot sites, checks their HTML, and exercises publication behavior.
It uses the existing licensed-check runner with `SYMBOLICA_LICENSE_SIGNED`; its
reusable Cargo artifacts are scheduled separately. `just ci-checks` selects the
same documentation check; uploads select its launcher and compiled producers.
Packaging, standalone Rustdoc, the persistent Typst renderer and WASM remain
additional checks in `nix flake check --impure`.

Compilation stays per crate. Source filtering preserves unaffected packages;
manifest and lockfile changes still invalidate broadly. Test dependencies avoid
multi-crate cycles: Linnet enables its own optional test features, and Spenso owns
the macro integration test. Regenerate Hakari after dependency/feature changes so
the shared cache does not retain unused build dependencies. Native package builds,
test harnesses and documentation share the optimized CI profile, compatible
features, Python interpreter and compiled library artifacts. Python extensions
retain their ABI and extension settings. Documentation stub exporters remain
separate because their global type registries can change the generated APIs when
combined. Unit-test crates and Rustdoc snippets still require their own compilation.
Clippy and doctests share check dependencies but retain workspace-wide source
inputs, so their rebuild isolation differs from package test artifacts.

Core API tests import UFO models using Python packages pinned in `uv.lock`.
These dependencies enter the runtime environment through `PYTHONPATH`; they do
not change the shared compiler interpreter or compiled test archives. The locked
Symbolica Python wheels support Linux on x86-64 and macOS on ARM. Requesting this
runtime on Linux ARM reports the missing upstream wheel.

The four NLO acceptance cases share one executable. Their assertions, workloads
and runtime worker reservations are unchanged.

== Fast tests and the full main suite

The shared Nextest configuration owns test selection. Ordinary local test
recipes and NixCI use `ci_gammaloop`: each test has a four-minute execution
limit, and a timeout fails. Compilation is outside that limit. Tests are
classified by measured duration, rather than excluding every `slow` module.
The measured extended scalar matrix contributes 139 cases to this fast lane;
the GL16 quartic case belongs to the long lane. A fast suite can still take
longer than four minutes because the limit applies to each test.

Every main push runs the supported full suite through `hydraJobsMain`.
Non-integration package groups use `ci_full`. Integration runs its ordinary
fast selection plus four disjoint `ci_slow` partitions, including all six
mandatory soft-counterterm acceptance cases. These partitions reuse the
existing optimized archives and do not introduce another feature set or
compile each shard independently. Missing main-run license credentials fail
the workflow rather than silently dropping licensed tests. Standard branch
Nix jobs continue to use `hydraJobs` and do not select the long partitions.

`ci_full` includes nonignored slow cases. Existing quarantines and manual
external-backend checks remain explicit exclusions in the configuration;
their absence from automatic coverage does not establish that they are
unsupported. Ignored tests are not enabled wholesale.

Each successful Nix runtime result retains per-package JUnit reports with
individual case durations. GitHub Actions exports test names, statuses and
durations to a CSV artifact. Captured test output, expressions and raw logs
are excluded. A failed Nix derivation has no registered report output; the
build still fails, even when there is no timing CSV to upload.
Nextest reservations continue to protect shared
generated-state folders and large evaluator allocations. Four shards permit
separate workers to run long cases concurrently; they do not remove these
protections within a worker. Partitioning by hash is deterministic but does
not guarantee an optimal balance of measured durations.

The Spenso group also runs `spynso3` unit and integration tests for typed tensor
APIs and shared Symbolica expressions. Spynso's embedded `typst/*.typ` files are
production inputs to its package and dependent builds. Its display integration
test mocks Typst and checks the missing-optional-renderer behavior.

Merged artifacts are self-contained compressed archives. Recursive inheritance,
writable extraction, Cargo fingerprints and epoch-1 timestamps must survive
changes to archive handling. Stable publication outputs depend on the artifacts,
not the commit. Lightweight test runners request the corresponding Nix check
result before fetching compiler state or test binaries. NixCI launches a cached
runner too, but a matching cached check result can avoid executing its tests.
The local upload command selects static checks, runtime launchers and compiled
producers, not runtime check results. Uncached runtime checks run in NixCI using
those prepared inputs. An empty successful test output alone does not retain
the necessary producers. Test groups are scheduled and reported independently.

The local upload command keeps all selected outputs rooted through temporary
result links until cache publication finishes. The command removes those links
on exit. This prevents concurrent garbage collection from deleting a realized
producer between the build and upload stages.

Synchronous dependency discovery incorporates Syd's PR \#104 suggestion while
retaining the explicit graph. Local uploads adapt \#105 through the Just command;
entering a shell installs no global upload hook.

== Performance measurements

`.github/scripts/run_performance.py` compares two normal CLI packages on one
worker, outside cached Nix check derivations. Build each revision once; every
measurement generates a fresh state and uses that revision's own serialized
format. The baseline and candidate commits, fixture, point corpus and settings
are recorded by their hashes. Both revisions receive identical settings and
run with one worker. Model and restriction files are imported from the shared
candidate repository and included in the workload hash; a builtin model name
would instead select each binary's own embedded model.

The Actions workflow authenticates NixCI downloads with the repository's
`NIXCI_NETRC` secret, containing only the `cache.nix-ci.com` machine entry.
It checks access before building, configures the Nix daemon to use that file,
and uploads both exact normal CLI closures after realization. Hestia remains
a read-only secondary cache. The temporary credential file has mode 0600 and
is removed on success or failure. Compiler products are cached; measurements
always run afresh.

The single manifest in `benchmarks/performance/suite.toml` defines four workloads
and three scalar routes, producing ten comparisons. Scalar workloads share their
point corpus and settings. The runner renders `scalar.toml` with each workload's
graph and route settings; `gg-hhh-gl15.toml` defines the physical workload.
The exact rendered card and declared inputs, including scalar DOT graphs and
shared model and restriction files, are hashed before measurement. Both
revisions receive those same bytes.

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Input], [Coverage], [Routes]),
  [Physical gg → hhh GL15], [Threshold-subtracted amplitude], [One parametric route],
  [Scalar GL04, q₁²], [Quadratic Lorentz numerator], [All three LU routes],
  [Scalar GL16, q₇²], [A different quadratic topology], [All three LU routes],
  [Scalar GL02, (q₁⁰)⁴], [Quartic energy numerator], [All three LU routes],
)

The scalar routes are orientation-local 3D with the summed evaluator,
explicit-sum 3D and projected local 4D with the single parametric evaluator.
All use local and integrated UV counterterms together with threshold
subtraction. The small DOT inputs preserve the source graph's edge ownership,
LMB, phases and vertex factors. Only the selected edge receives the numerator
probe and its corresponding degree of divergence. Both revisions import the
same candidate-side graph bytes, included in the fixture fingerprint; UV
forests, counterterms and evaluators are generated afresh in every trial.
These fixtures exercise numerator and route differences without relying on
the exactly vanishing massless cuts of several Lorentz-quartic examples.

Generation measures the complete import, generation and state-save command.
For scalar cases, graph enumeration is outside this timing because the graph
is supplied as a fixed DOT input. GL15 also runs graph enumeration.
Numeric metadata retains the native per-graph total and the Spenso, Symbolica
and compilation durations from `generation_summary.json`. These explain
generation costs without changing the wall-time metric or its thresholds;
the native total excludes CLI initialization, import and state saving.
Runtime measures the existing `bench` command at four distinct fixed points,
with caching disabled and evaluator compilation disabled. Each point is
repeated within its own warmed benchmark. GL15 fixes one helicity/external
configuration and selects one orientation and LMB channel, with stability
rotations disabled. Scalars retain the correctness tests' rotation settings
and sum all orientations and LMB channels at nine-dimensional points.
`bench` disables event generation for its timed batches; selectors and
observables are absent. This measures fixed-point integrand evaluation rather
than a batch of changing integration points. Large scalar numerators,
compiled evaluators and isolated cold-load cost remain outside this fast gate.
The extended scalar correctness suite supplies timing telemetry for the
long Lorentz-quartic cases, including GL16 q₇⁴, on main.
The requested benchmark duration per point is five seconds for GL15 and two
seconds for scalars. The CLI calibrates sample counts from its ten-sample
warmup, so the actual warmed loop can be shorter than that target.

Five baseline/candidate pairs alternate execution order. A suspected regression
requires another five pairs starting in the opposite order. Generation must
increase by more than 20% and 0.5 seconds; warmed runtime must increase by more
than 15% and one microsecond per sample. At least four pairs must be slower,
and the effect must exceed three times the robust paired noise estimate.
This is a noise-screening policy, not a statistical confidence interval.
Every physical result must be finite, nonzero and agree between revisions.
Scalar results use relative tolerance 10⁻⁹ with no absolute allowance, because
an absolute floor suitable for GL15 could mask changes to their smaller values.
Each topology/route comparison, including any confirmation run, has a
four-minute budget. The collection can take longer than four minutes.
Compilation is outside that budget. Exceeding it yields
an inconclusive failure rather than a claimed performance regression.

The comparator exits zero for a clean comparison, one for a confirmed
regression or changed result, and two for invalid, incomplete or inconclusive
measurements. The public output directory contains revision/workload
fingerprints, numeric measurements and the comparison report. Generated states
and raw process logs stay in the private work directory.

Uploading these measurements does not enforce a performance budget. A merge
gate must run the harness freshly and require its zero exit status. The existing
required `NixCI readiness` check aggregates the performance job together with
configuration validation. A planned measurement must succeed; failed,
cancelled or unexpectedly skipped measurements fail the aggregate.

Automatic paired benchmarks run for non-draft, same-repository PRs marked
`final-review` when the pinned candidate/base diff changes Rust code, model
or calculation inputs, dependencies, build settings or the benchmark itself.
Documentation-only and integration-test-body changes skip the measurement.
Merge-group commits apply the same relevance test using the queue's head/base
SHAs. Main pushes explicitly skip paired benchmarks: their full correctness
suite and timing uploads continue independently. Ordinary feature-branch pushes
do not trigger this workflow. Manual benchmark dispatch remains available.

A relevant new commit while `final-review` remains set must rerun the gate.
Returning to development means removing that label; PR-number concurrency
cancels superseded work. Compiled CLI packages are reused from Nix caches, but
the benchmark always measures fresh states. All ten cases run sequentially on
one worker in a single runner invocation and pinned runtime shell, to avoid timing
interference and reuse the two CLI packages built once before measurement.
The runner gives each case separate work and output subdirectories. Every case
must pass; a failing case does not prevent the remaining cases from recording
their own reports. The performance manifest and cards are excluded from the
drawing WASM source, which the CLI embeds, so changing a benchmark does not
invalidate that compiled dependency. Runtime tools come from the pinned Nixpkgs
input without building the full development shell.
The signed license is scoped to the measurement step, and uploads contain
only numeric evidence and revision/workload metadata.

== Final-review readiness

`nix-ci.nix` is self-contained. Its top-level `enable` boolean is manually
controlled; `just ci-update` reads and preserves that value while regenerating
all other fields. The graph check still compares the complete generated file,
so either boolean is valid locally but missing/nonboolean values and scheduling
drift fail. `dependency-discovery.enable` remains independent. Disabling remote
NixCI does not disable local checks or cache uploads.

The NixCI readiness workflow first evaluates the PR head configuration.
Its `NixCI readiness` check requires a non-draft PR with `final-review`
and top-level `enable = true`. Draft and unlabeled PRs fail this merge gate while
remaining usable for development and feedback. Main pushes and merge-group
commits require enabled configuration without a label condition.

The full GitHub Actions build workflows run automatically on pushes to `main`,
not on PRs, merge groups, or feature-branch pushes. This includes Nix, Continuous integration,
Documentation Pages, Linnet Python WebAssembly, and Typst package mirroring.
Existing manual dispatch, scheduled maintenance/acceptance checks, and release-tag
publishing remain available. The readiness workflow retains PR and merge-group
triggers and invokes the scoped performance job described above. The
`final-review` label gates readiness and PR performance measurements, not the
full build workflows or NixCI itself; labels
do not alter committed NixCI configuration or cancel already-running NixCI jobs.

Agents disable NixCI during implementation unless instructed otherwise. Before
final review, enable it, validate/upload the final code, push, and then apply the
label. Returning to development means removing the label and committing
`enable = false`. Keep `main` enabled. See
#link("../../CONTRIBUTING.typ#ci-readiness")[the contributor workflow].

For rollout, create the `final-review` label first. After the readiness workflow
is merged and its successful check is registered, add `NixCI readiness`, bound
to GitHub Actions, to the existing `PRs must pass` ruleset. Preserve the required
`deploy packages.x86_64-linux.nix-ci-passed` status from NixCI and all unrelated
rules. The readiness check supplements successful tests; it does not replace them.

Prefer an available completion notification over keeping an agent polling CI.
Retain the commit/run identity and logs, yield after dispatch, and verify the
result when notified. Local validation and uploads must finish successfully
before pushing. Notification support depends on the execution environment;
this repository does not install a CI-to-agent wake-up service. See
#link("../../CONTRIBUTING.typ#ci-completion")[the completion workflow].

== Branch maintenance

After rebasing an active branch onto the CI changes, preserve its test groups and
source filters, run `just ci-update`, and commit both generated files. FeynKit
needs its eighth group, nine extracted crates, embedded fixtures, optional Python
dependencies, boundary checks and Python stub checks. Its extracted libraries
must remain independent of GammaLoop and its private workspace-hack dependency.
The shared FeynKit branch was not changed by the experiments.

The `ci-cache-base` input pins a compatible compiler-state seed. Keep that commit
reachable in published history. Refresh deliberately after a green revision with
`just ci-cache-base FULL_COMMIT_SHA`; rebasing or cleaning history alone is not a
reason to invalidate it. See the #link("nix-crane-cache-reuse.typ")[cache-reuse audit]
for earlier investigation and implementation details.

== Measured results

The Cargo dependency cleanup was compared locally with `246c03693` on Rust
1.98.1, using eight jobs, `ci-optim`, no incremental compilation, and the same 631
selected tests from the Linnet and Spenso groups (658 listed, identical filters).

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Local Cargo workflow], [Before], [After]),
  [Empty-target test build and execution], [6m52s], [6m13s],
  [Empty-target check followed by tests], [10m31s], [7m12s],
)

These are single pairs on a shared host with vendored dependencies available,
not full-suite or remote measurements. Direct test CPU time fell 4.4%; most of
the larger combined gain comes from checking. Both unchanged retries compiled
nothing. Nix archive inventories matched all 631 tests; the selected Nix suite
passed 1,633 runtime tests (including five Python cases), 43 doctests and static
checks. FeynKit graph evaluation preserved its groups and extracted-crate cache
boundaries; its runtime suite was not rerun for this cleanup.

The final warm comparison used the same main application layout:
#link("https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-clan-cache-fixed-validation/f1ed70f585d293ecb851cf39f28821a69b33ac05")[baseline]
and #link("https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-lazy-check-runners/dd6ba768692302ca4855727110c4d8422b0979f9")[candidate].

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Warm NixCI measurement], [Before], [After]),
  [Required checks], [8m26s], [3m46s],
  [Whole suite], [8m37s], [4m05s],
  [Observed worker minutes], [6.98], [3.96],
  [Logged artifact restores], [1.31 GiB], [8,160 bytes],
  [Final success after required checks], [11s], [19s],
  [Cached build jobs], [53/53], [53/53],
)

All eight test jobs passed without compilation or executed test summaries. Their
truncated output paths prevent the generic reporter from proving result reuse
for every job. This is one warm pair; waiting conditions also differed. Reported
artifact sizes exclude other traffic, including 13.4 MiB shown by one inner Nix
command. Worker minutes are observed log spans, not billed time or CHF.

Two unresolved observations accompany this pair: 46 previously readable unchanged
cache roots returned HTTP 404 before republishing, and the doctest appeared started
with empty logs before a 78-second gap in observed worker activity ended. Their
causes are unknown. No abandoned workers, repeated builds, replaced logs or
transfer errors were observed in the completed pair.

Earlier local comparisons measured different scopes:

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Controlled local runtime comparison], [Before / Cargo], [After / Nix]),
  [Main source edit, including evaluation/setup], [10m12s], [4m53s],
  [FeynKit source edit, including evaluation/setup], [9m58s], [5m00s],
  [Main project-cold build and tests: Cargo versus Nix], [21m46s], [28m48s],
)

These include Python preparation and all runtime groups, but exclude Clippy and
doctests. Each comparison matches application code within its own layout; it does
not measure the crate split itself. The cold pair used identical toolchains and
eight CPUs, an empty Cargo target and a disposable Nix store without project
outputs. External toolchains, native libraries and crate sources were available;
remote builders, downloads and uploads were disabled. Nix still took 32.3% longer.

The complete local correctness runs passed 1,633 tests/43 doctests on main and
1,834 tests/62 doctests on FeynKit, with matching inventories, filters and existing
skips, including five Python cases per layout. FeynKit's boundary checks passed;
its Python-stub workflow was preserved but stub generation was not rerun. A
pre-existing probe expecting a `feynkit-py` edit to affect only that crate failed:
its API consumer also depends on it. That assertion was not changed. An unchanged
retry passed two intermittent main state-file failures; that failed attempt is
excluded from successful timings. The final runner fix separately passed eight
cached checks without compiler artifacts/downloads, a real five-test Clinnet cache
miss, and a missing-license failure. Consolidation preserved all 63 compared Nix
output identities and all 11 required local checks; those checks reused results.

Detailed evidence remains in these file histories, ending at the pre-cleanup
commit `0cc8cd2`:
#link("https://github.com/alphal00p/gammaloop/commits/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-final-warm-comparison.json")[final warm data history],
#link("https://github.com/alphal00p/gammaloop/commits/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-artifact-preparation-results.json")[local artifact and coverage history],
#link("https://github.com/alphal00p/gammaloop/commits/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-cold-cargo-comparison.json")[cold commands and results history], and
#link("https://github.com/alphal00p/gammaloop/commits/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs")[experiment history].
This keeps generated reports and experiment diaries out of the maintained source tree.

== Collecting another comparison

Run `just ci-report MANIFEST OUTPUT_DIR`. It selects pinned Node and uses an
authenticated `gh` for GitHub. The reporter only reads existing runs; it never
launches, retries or cancels jobs. Put generated output outside the source tree.
A minimal manifest for reading an existing suite is:

```json
{
  "repository": "alphal00p/gammaloop",
  "runs": [{
    "layout": "main",
    "variant": "candidate",
    "scenario": "warm",
    "sha": "dd6ba768692302ca4855727110c4d8422b0979f9",
    "suiteUrl": "https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-lazy-check-runners/dd6ba768692302ca4855727110c4d8422b0979f9"
  }]
}
```

Add a `baseline` entry with matching `layout`, `scenario` and optional `pair`
(default `"1"`) for paired percentages. Use a new pair label for repeats. Each
entry requires a full SHA and exact commit-level suite URL. Duplicate suites and
ambiguous pairs are rejected. Optional `actionRunIds` limits workflows; otherwise
all Actions runs for the SHA and every attempt are retained.

For controlled comparisons, supply the same explicit `requiredAttributes` in both
entries: every test group, doctest, Clippy, formatting and graph check. Otherwise
the default is all observed test jobs plus the three static checks. A runner's
build alone does not satisfy its required test phase. Missing groups, incomplete
logs, failed suites or ambiguous clocks prevent accepted paired percentages.
Coverage equivalence is never inferred: compare inventories, filters, features,
Python behavior, license settings and static-check results separately.

The output contains `report.md`, `report.json`, CSV tables for suites, jobs,
Actions, transfers, compilations, repeated derivations, incidents and pairs, plus
`raw/N/` evidence and the input manifest. Use a fresh directory for each snapshot;
existing files are overwritten. New raw files are owner-readable/writable; review
them before sharing. Exit status is 0 for complete evidence, 2 for a partial report,
and 1 for invalid invocation or manifest.

Metric interpretation:

- Required latency ends at the last intended check; whole-suite latency ends at
  the last completed job. Final-success tail is reported separately. Recorded
  `submittedAt`/`submissionCompletedAt` push bounds can additionally measure elapsed
  time from submission; commit author dates are not submission clocks.
- Worker minutes sum observed worker spans. Check duration and pre/post-worker
  gaps stay separate; a pre-worker gap does not prove queue time.
- Only completed transfer messages supply sizes/times. Intervals are unioned per
  worker before summing, because downloads and uploads can overlap. Content sizes
  are not wire bytes; timers may include decompression, import and locking.
- `Compiling`, `Checking` and Cargo command durations remain separate by context.
  Repeated crate names can mean different features or profiles. Only the same
  full derivation path built under different job URLs establishes repeated Nix work.
- Test execution is `executed`, `reused` or `unknown`. Fetching a runner wrapper
  does not prove reuse of its test result. Nextest `failed` includes `timedOut`;
  do not add those counters together.

To replay without network calls, add `offline` to an entry, with `suite`, `checks`,
`jobs` and `actions` paths pointing at a saved `raw/N/` snapshot's corresponding
JSON files. Paths resolve relative to the manifest; log paths in `jobs.json`
resolve relative to that index. The reporter accepts the earlier collector's job
index too. Keep raw evidence alongside the report so replays remain possible.

To preserve unusual service behavior, retain snapshots when observed; final logs
may overwrite earlier attempts. The optional manifest fields are:

#table(
  columns: (1fr, 1fr),
  table.header([Field], [Evidence to save]),
  [`observationFiles`], [Snapshot JSON with `at` and `suites: [{spec: {sha, suiteUrl}, suite: API_RESPONSE}]`],
  [`dependencySnapshot`], [`{sha, file}` for that exact revision's generated config containing `dependencies`],
  [`observedInterruptions`], [Entries `{jobUrl, observedAt, evidenceFile}`; evidence repeats the URL with `status: "abandoned"` and optional prior log `file`],
  [`observedLogReplacements`], [Same entry fields; evidence repeats URL/time with earlier `file` and later `replacement.file`],
)

Evidence paths resolve relative to the manifest, linked logs relative to their
evidence file. The reporter copies them into `raw/N/` for offline replay; point
future manifests at those copies. Replaced logs and abandoned work make current
resource totals lower bounds, even after a successful retry. Prior totals are not
added because streams may overlap. Queue observations require every declared
prerequisite to be successful in the same snapshot; they do not prove continuous
queue time or its cause. Transfer errors, invalid clocks, mismatched outcomes and
repeated attempts retain evidence without attributing a backend cause.

Run the reporter's offline regression checks with
`node --test .github/scripts/ci-report.test.mjs`.

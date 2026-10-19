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

Actions use Rust 1.98.1, matching the Nix toolchain, for native and WebAssembly
builds. Update the toolchains deliberately together.
Main's native Actions jobs install Just 1.58.0 rather than the Ubuntu package,
which cannot parse the shared recipe attributes. Fast-test jobs also install
Typst 0.15.1 and prebuilt Jujutsu 0.45.1 for drawing and provenance tests.
These tools are verified before the test suite starts on Linux and macOS.
Manual dispatches can set `run_native_tests = false` for static-check debugging.
It defaults to true; main pushes always run the native tests. Nix's runtime suite
is independent of this manual option.

Nix Actions producers and matrix workers authenticate to the NixCI cache before
realizing their exact inputs. They build one derivation at a time with two build
cores and retain output roots through Hestia publication. Nix's automatic store
collection keeps a 4–8 GiB free-space watermark; the producer checks for at least
12 GiB free before compiling and reports free space after realization. The
8 GiB producer-publication limit still applies to the final closure, independently
of temporary compiler/store usage.

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

Prebuild and normal dependency commands omit feature selections already encoded
by their generated Cargo anchors. Test dependency commands retain features on
selected component crates and all forwarded selectors, since a selector can
activate an optional dependency on a selected root. Expanded feature vectors still define
cache contexts and workspace checks.

On platforms with a locked Symbolica wheel, core API tests import UFO models
using Python packages pinned in `uv.lock`. These dependencies enter the runtime
environment through `PYTHONPATH`; they do not change the shared compiler
interpreter or compiled test archives. The wheel tags in `nix/ci.nix` define the
supported platforms: Linux on x86-64 and macOS on ARM. On Linux ARM, the core
group omits this runtime and excludes the two UFO-import tests through a Nextest
run-time filter, while loader-version validation still runs. Archive filtersets
cannot select individual tests. Explicitly requesting
`gammaloop-ufo-python-runtime` on Linux ARM reports the missing upstream wheel.

The NLO acceptance cases share one executable. Their assertions, workloads
and runtime worker reservations are unchanged.

== Fast tests and the full main suite

The shared Nextest configuration owns test selection. Ordinary local test
recipes and NixCI use `ci_gammaloop`: each test has a four-minute execution
limit, and a timeout fails. Compilation is outside that limit. Tests are
classified by measured duration, rather than excluding every `slow` module.
Other scalar cross-section tests remain eligible despite their `slow` module
names; the GL16 quartic case belongs to the long lane. A fast suite can still take
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

GammaLoop's adapter in `.github/scripts/run_performance.py` prepares one suite per
manifest case and delegates paired execution, numerical qualification, regression
screening and sealed evidence to the commit-pinned gammabench framework.
The fixture manifest and native measurement semantics remain in this repository.
`benchmarks/performance/gammabench.rev` pins the framework by complete commit.

`.github/scripts/run_performance.py` compares two normal CLI packages on one
worker, outside cached Nix check derivations. Build each revision once; every
measurement generates a fresh state and uses that revision's own serialized
format. The baseline and candidate commits, fixture, point corpus and settings
are recorded by their hashes. Both revisions receive identical settings and
run with one worker. Scalar model and restriction files are imported from the
shared candidate repository and included in the workload hash. The SM cases
use the builtin `sm-default` name instead, so each binary loads its own
embedded model and a change to the SM JSON format cannot stop the baseline
from parsing it. Every result must still agree between revisions.

The Actions workflow authenticates NixCI downloads with the repository's
`NIXCI_NETRC` secret, containing only the `cache.nix-ci.com` machine entry.
It checks access before building, configures the Nix daemon to use that file,
and uploads both exact normal CLI closures after realization. Hestia remains
a read-only secondary cache. The temporary credential file has mode 0600 and
is removed on success or failure. Compiler products are cached; measurements
always run afresh.

The single manifest in `benchmarks/performance/suite.toml` defines thirteen
CFF revision comparisons: the physical workload, nine scalar route cases,
two evaluator-layout cases and one thermal case. A fourteenth fast comparison
runs LTD and CFF on the same candidate; the slow Arb comparison requires
explicit selection. Scalar route cases share their point corpus and settings.
The runner renders `scalar.toml` with each workload's graph and route settings;
the other cases use fixed cards. The exact rendered card and declared inputs,
including DOT graphs and shared model and restriction files, are hashed before
measurement. Both
revisions receive those same bytes.

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Input], [Coverage], [Routes]),
  [Physical gg → hhh GL15], [Threshold-subtracted amplitude], [One parametric route],
  [Scalar GL04, q₁²], [Quadratic Lorentz numerator], [All three LU routes],
  [Scalar GL16, q₇²], [A different quadratic topology], [All three LU routes],
  [Scalar GL02, (q₁⁰)⁴], [Quartic energy numerator], [All three LU routes],
  [Scalar GL16 evaluator probe], [Reachable aliases and selected function maps], [All four evaluator methods],
  [Two-loop gluon self-energy GL1], [Weighted local-UV guards and runtime arithmetic], [All four evaluator methods],
  [Thermal gluon self-energy GL1], [QCD numerators with thermal weights and nested local UV forests], [SingleParametric and SummedFunctionMap],
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

The two evaluator probes exercise complementary performance mechanisms.
GL16 includes local and integrated UV terms, with 174 orientations in each of
two evaluator stacks and substantial shared alias catalogs. GL1 retains nine
orientations with one fermion-edge sign fixed, plus local UV terms. Integrated
addbacks, soft counterterms, threshold subtraction and thermal settings are
disabled for GL1. Separate retained-source validation found repeated weighted
guards in GL1: shared regions reduced root guards from 93 to 21 and final
SingleParametric root inversions from 639 to 299. These diagnostic exports are
outside the measured path. Both probes build SingleParametric, Iterative,
SummedFunctionMap and Summed, then check all four at the same two fixed points.
They pin graph and LMB selectors and sum their retained orientations, with
stability rotations disabled.
The additional `gl04-q1-squared-ltd-vs-cff-4d` case uses the same scalar card,
four points, local 4D UV counterterms and threshold subtraction. It generates
CFF and LTD separately with the same candidate executable, and selects each
representation explicitly at a single Double stability level. Every timed
sample must report exactly one successful two-probe attempt in that representation.
After each timed batch, displayed real and imaginary extrema are collected using
the same per-sample Jacobian normalization as the warmup. The LTD/CFF case rejects
invalid values or extrema outside the existing component-wise correctness
tolerance; this check adds no work to the measured interval.
Native residue-map sizes are recorded separately; physical cut and threshold
edge identities must match. This is a same-revision comparison, not a historical
LTD measurement on a baseline that does not support LTD. Its two rendered-card
hashes differ, while the shared physical settings, corpus, binary and commit
identities must agree. Both generation and warmed runtime use a 50% relative
ceiling, retaining the 0.5-second and one-microsecond absolute floors and the
same paired-noise and independent-confirmation rules. The thirteen CFF revision
comparisons keep their original policies.

Only this frozen GL04 q1² comparison has an imaginary-zero certificate. The
unique simple massless edge 1 is canceled by its numerator; both cuts carrying
threshold counterterms contain that edge. The surviving cut `[3,4]` has no
thresholds and a real scalar phase, preserved by the MUV/MSbar integrated UV
normalization and real sampling/flux factors. The producer checks frozen source,
model, card and corpus identities, saved phase-relevant defaults, and the actual
simple-cut/threshold catalogue. Each of the four imaginary remainders must
independently satisfy `abs(Im) ≤ sqrt(2⁻⁵²) E_cm³` at the fixed `E_cm = 1`, using
the existing dimensioned scalar zero policy. Real components retain the 10⁻⁹
relative check with no absolute allowance. Raw outputs and strict mismatch counts
remain recorded; timed extrema retain their strict per-point warmup comparison.

The thermal probe generates the same GL1 graph, without its orientation
pattern, in thermodynamic equilibrium at inverse temperature one and muB = 3,
with a fixed timelike external momentum. Its quark loop and triple-gluon vertex
supply QCD numerators with thermal weights. Vacuum subtraction is off: it
would remove the two-loop forest in which both overlapping one-loop
counterterms nest. Only local UV counterterms are built; integrated addbacks
more than triple its generation time, and building all four evaluators would
also exceed the case budget. It builds and checks SingleParametric and
SummedFunctionMap at two fixed points, with graph and LMB selectors pinned,
all orientations summed and stability rotations disabled. Three-loop QCD
vacuum graphs keep nested forests with vacuum subtraction, but each takes
minutes and over 1 GB to generate on one worker.

Generation uses Hyperfine 2.0.0 to measure the complete import, generation and state-save CLI process. Its raw export is retained and verified against the numeric metric; this scope supersedes the earlier Python-controller clock.
For scalar cases, graph enumeration is outside this timing because the graph
is supplied as a fixed DOT input. GL15 also runs graph enumeration.
Numeric metadata retains the native per-graph total and the Spenso, Symbolica
and compilation durations from `generation_summary.json`. The native total
excludes CLI initialization, import and state saving. Probe cases also gate
each evaluator's construction duration from the existing structured milestones,
summed across the graph's complete evaluator stacks. Missing, duplicate or
unfinished stages invalidate the measurement. Whole-generation peak RAM is
gated separately; it is not attributed to individual evaluator modes.
Runtime measures the existing `bench` command at four distinct fixed points
for the original cases, and two points for each evaluator or thermal probe,
with caching disabled and evaluator compilation disabled. Each point is
repeated within its own warmed benchmark. GL15 fixes one helicity/external
configuration and selects one orientation and LMB channel, with stability
rotations disabled. Scalar route cases retain the correctness tests' rotation settings
and sum all orientations and LMB channels at nine-dimensional points.
`bench` disables event generation for its timed batches; selectors and
observables are absent. It resolves the model and warms up the integrand once
before its warmup samples, so neither the calibration nor any timed batch
includes that preparation. This measures fixed-point integrand evaluation rather
than a batch of changing integration points. Large scalar numerators,
compiled evaluators and isolated cold-load cost remain outside this fast gate.
The extended scalar correctness suite supplies timing telemetry for the
long Lorentz-quartic cases, including GL16 q₇⁴, on main.
The requested benchmark duration per point is five seconds for GL15, two
seconds for the scalar route cases, 0.75 seconds per evaluator for GL16 and
thermal GL1, and one second per evaluator for GL1. Shorter runs let runner
timing noise alone make unchanged warm runtimes inconclusive, while a
comparison with confirmation must still fit its four-minute budget. Every point uses ten warmup samples and ten timed batches. The CLI
calibrates sample counts from its warmup, so the actual warmed loop can be
shorter than that target. Each revision generates its state once per trial,
regardless of the number of methods.

Measurement schema 2 records a separate numeric output vector and runtime
metric for every requested evaluator method. Each method must match its own
revision's reference evaluator and the baseline reference at every corpus point.
This also rejects a method-specific incorrect result shared by both revisions.
The runner generates one fresh state per revision and trial, then switches
evaluator methods in one read-only benchmark session. Timings are never
averaged across evaluator methods.

Five baseline/candidate pairs alternate execution order. A suspected regression
requires another five pairs starting in the opposite order. Generation must
increase by more than 20% and 0.5 seconds; each evaluator's construction by
more than 20% and 0.1 seconds; each warmed runtime by more than 15% and one
microsecond per sample; peak RAM by more than 20% and 64 MiB. At least four
pairs must be slower, and the effect must exceed three times the robust paired noise estimate.
Noise can make a comparison inconclusive only when the paired absolute delta
plus its noise margin also exceeds the absolute budget. Large relative noise
on a tiny operation does not by itself block the gate.
This is a noise-screening policy, not a statistical confidence interval.
Every physical result must be finite, nonzero and agree between the paired
revisions or representations.
Scalar results use relative tolerance 10⁻⁹ with no absolute allowance, because
an absolute floor suitable for GL15 could mask changes to their smaller values.
Each topology/route comparison, including any confirmation run, has a
four-minute budget. The collection can take longer than four minutes.
Building the two CLI binaries is outside that budget. Exceeding it yields
an inconclusive failure rather than a claimed performance regression.

The explicit `paper-massless-soft-cff4d-arb` case uses the existing paper B.1
massless soft-counterterm graph and its complete CFF/from4D sum. It retains
the original IR prescriptions, external momenta and masses; that fixture
disables integrated UV and threshold counterterms. The stored binary64
momentum bits come from the original full-subset UV-profile ray with seed
1342 and scale 10¹². Complete signed edge routing certifies each revision's
positive-unit permutation into the canonical physical loop-edge basis.
Catalogue ordering, internal tree serialization and historical timing values
are not part of this certificate.

Numerator inlining is explicitly enabled, matching the production default.
The benchmark sets `store_atom = false` to avoid retaining symbolic expressions
that are unnecessary for numerical evaluation.

This slow gate generates one fresh state per revision, then runs five AB/BA
runtime pairs through the same shared `precise_sample.rs` producer, compiled
independently against each revision's normal API through its existing Crane
crate-build owner. Each trial starts a fresh process, resolves the integrand's
model, calls public runtime warmup, performs one numerical warmup and measures
one forced 1000-bit Arb call with point caching disabled. Both calls must have
one stable Arb attempt with the original two probes. The gate checks both native
`evaluator_evaluation_time` and the complete warmed physical call, so moving work
outside the evaluator timer cannot hide a regression. State loading and model
warmup are recorded separately. Generation is a single preparation observation,
not a paired generation verdict. The thirteen CFF revision cases' generation policy is
unchanged. A declared two-hour budget includes preparation and any independent
timing confirmation. The same 15% plus one-microsecond runtime policy applies.

The slow gate preserves both native decimal components and the strict 10⁻⁵
relative comparison as diagnostic evidence. The fixture-specific source predicate
certifies only the real component as analytically zero;
when both producers pass that predicate, each real remainder must independently
satisfy `abs(Re) ≤ sqrt(2⁻⁹⁹⁹) E_cm⁻⁷`, following the existing dimensioned scalar
zero policy. The imaginary comparison retains the original relative tolerance
with no absolute allowance. No binary64 conversion or other component's scale
enters these comparisons. Raw strict mismatch counts remain in the report.
A numerical mismatch remains a failure while timing measurements and required
confirmation continue; a timeout remains incomplete evidence. Only this precise
case sets `confirm_on_correctness_failure = true`; the fast cases retain the
framework's default confirmation behavior.

These two certified-zero cases use gammabench's explicit `certified_zero`
contract: a certificate with an `id` and `sha256`, and component-specific
`index` and `absolute_bound` entries. Each adapter result must carry the exact
declared certificate. GammaLoop's `GLFixture` owns the physical proof and
dimensioned bounds; the framework verifies certificate identity and applies
only the declared component allowances. It retains raw values and strict
mismatch diagnostics alongside the qualified result. The slow case selects
`correctness_format = "decimal"` so this comparison never rounds through
binary64, and `budget_class = "slow"` for its 7200-second deadline. Ordinary
cases retain the 240-second default. Each role prepares its state once inside
its first trial, sharing that deadline with all paired samples and any
confirmation. This preparation is not a paired generation metric.

Use an adequately provisioned host such as `itphlies`: exact-main generation
has exceeded 21 GB, whereas the
#link("https://docs.github.com/en/actions/reference/runners/github-hosted-runners")[standard hosted Ubuntu runner]
provides 16 GB.
The existing Actions job explicitly selects `--lane fast`; it does not claim
coverage of the slow gate. No larger paid or self-hosted runner is configured
by this change. To run the slow gate locally, build the normal CLI packages and
the two shared precise producers at the actual baseline and candidate commits:

```sh
ulimit -s 65536
export RUST_MIN_STACK=67108864
nix build --impure --file benchmarks/performance/precise_sample.nix \
  --argstr revisionRoot /path/to/baseline --argstr revision "$BASELINE" \
  --out-link /tmp/precise-baseline
nix build --impure --file benchmarks/performance/precise_sample.nix \
  --argstr revisionRoot /path/to/candidate --argstr revision "$CANDIDATE" \
  --out-link /tmp/precise-candidate
python3 .github/scripts/run_performance.py \
  --case paper-massless-soft-cff4d-arb --input-root "$PWD" \
  --baseline-commit "$BASELINE" --candidate-commit "$CANDIDATE" \
  --baseline-binary "$BASELINE_CLI" --candidate-binary "$CANDIDATE_CLI" \
  --baseline-precise-producer /tmp/precise-baseline/bin/precise_sample \
  --candidate-precise-producer /tmp/precise-candidate/bin/precise_sample \
  --work-dir /fresh/private/paper-arb --output-dir /fresh/public/paper-arb
```

Both producers receive identical source and fixture bytes. The manifest hashes
the physical certificate, graph, model, prescriptions and producer source;
measurement metadata also records the actual executable hashes. Run with the
usual licensed environment and memory watchdog. `--lane all` includes this
slow case together with the fourteen fast comparisons.

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
SHAs. The separate `Performance main` workflow measures every main push, including documentation-only merges, against the previous main SHA. These captures run independently of PR and merge-group gates. Ordinary feature-branch pushes
do not trigger this workflow. Manual benchmark dispatch remains available.

A relevant new commit while `final-review` remains set must rerun the gate.
Returning to development means removing that label; PR-number concurrency
cancels superseded work. Compiled CLI packages are reused from Nix caches, but
the benchmark always measures fresh states. All fourteen fast cases run sequentially on
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

=== Sealed results and premerge previews

Each case retains discovery and any confirmation as separate immutable snapshots.
The five-pair confirmation begins in the opposite order and shares the original
240-second case budget. Numeric evidence includes the complete source commits,
Nix package identities, fixture hashes, output vectors, native metrics and
Hyperfine exports. Private states and logs are never part of the upload.

`Publish performance results` runs from the default branch after the measurement
workflow completes. It verifies the source run, PR identity and artifact before
using `GAMMABENCH_PUBLISH_TOKEN`, a fine-grained token restricted to gammabench
Contents read/write. Candidate code is not executed in that publisher.
Main captures append to gammabench main; PR captures append to its separate
`pr-results` branch. One Pages deployment combines main history with
#link("https://alphal00p.github.io/gammabench/previews/")[PR previews].
Each PR has a stable landing page and permanent links for each measured SHA,
run and attempt. PR evidence never becomes a main-history point on merge.

Failed and inconclusive measurements are published too. Publication errors are
separate from the numerical gate. To retry publication without rerunning the
benchmark, manually dispatch `Publish performance results` with the original
Actions run ID and attempt while its 14-day artifact is available. Identical
redeliveries are no-ops; conflicting evidence is rejected. Results committed to
gammabench have no automatic expiration.

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
`NixCI` suite status from NixCI and all unrelated
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

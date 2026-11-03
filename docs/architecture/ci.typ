= CI maintenance and measurements

Use `just ci-checks` for the selected local suite and `just ci-checks-and-upload`
before final review. The latter runs static checks and publishes compiled runtime
producers without running the runtime suite. NixCI runs that suite after the push.
See #link("../../CONTRIBUTING.typ#nixci-cache")[CONTRIBUTING.typ] for credentials,
licenses and retries. `just check` checks compilation through the selected
backend without running tests.

== Native development commands

The root Just recipes default to `backend=auto`: native tests, builds and PyO3
packaging use Nix compilation and runtime runners when Nix is available,
otherwise Cargo/Nextest or Maturin. `scripts/develop.sh` resolves and prints the
backend before starting backend tooling. Without Nix, Rust commands invoke
Cargo and Nextest directly, with no Python dispatcher. With Nix it obtains
pinned shell tools from `flake.lock`, then `scripts/develop_nix.sh` invokes the
native build derivations and test runners. Python remains in the separate wheel
and cache-inspection tools. A Nix failure is never retried through Cargo.
Use `just backend=cargo COMMAND`, `just test-cargo CRATE` or
`just build-cargo CRATE` for explicit local Cargo iteration. An explicit Cargo
override uses host tools. Use `nix develop -c just backend=cargo COMMAND` when
local compilation needs the pinned toolchain and native libraries. See
#link("../../CONTRIBUTING.typ#development-commands")[Common Commands] for command
scope, specialized workflow exceptions and local `fmt` / `clean` operations.

`backend=dyndrv` is a separate explicit choice, using the `cargo-dyndrv` fork
pinned in `flake.lock` at `c1935d6`. It shares Cargo command construction with the
ordinary Cargo path and delegates source staging, compiler delivery and store
lifecycle to the pinned worker, with the local native runtime-path patch applied
by `nix/dyndrv.nix`. It does not replace the NixCI graph or upload
selection. See #link("../../CONTRIBUTING.typ#dyndrv-development")[the dyndrv
workflow] for host requirements, qualification commands and retained state.

The Nix shell backend snapshots paths tracked by JJ or Git and copies
their live contents, including uncommitted edits. It excludes `target/` and
`.venv/`, and ignored files are not imported just because they exist beneath
the workspace root. Newly added inputs must be tracked. The bootstrap pins
shell tools before constructing the snapshot, keeps its tooling and source alive
for the command, and cleans up afterward. Nix evaluation and
build subprocesses receive an environment without `SYMBOLICA_LICENSE` or
`SYMBOLICA_LICENSE_SIGNED`.

Use `just test CRATE [NEXTEST_ARGS...]` for an individual workspace crate,
including libraries outside GammaLoop's core. `just test-all` and `just test-ci`
retain whole-workspace selection. For example:

```sh
just test spenso
just test linnest -E 'test(layout)'
just test gammaloop-integration-tests --test test_runs
```

Native commands share `.config/build-profile` (default `ci-optim`, matching CI).
`just profile=dev COMMAND` and `just profile=release COMMAND` select explicit
profiles for Rust and Python alike. The environment override is
`GAMMALOOP_BUILD_PROFILE`. For tests, `--cargo-profile` changes compilation
and `--profile` selects the Nextest runtime profile.

Cargo uses one native Nextest invocation with `--workspace` or repeated `-p`
selections. Features follow Cargo defaults and the caller's explicit flags;
`--lib` and `--test` are Cargo build selectors. This avoids repeated compilation
across per-package feature contexts and produces one JUnit report containing the
selected tests, without later package runs overwriting earlier results. Cargo
listing emits native Nextest JSON directly. Python API tests remain opt-in:
`just test-python-api` adds `--features python-api-tests` and requires a prepared
Python environment, as described in
#link("../../CONTRIBUTING.typ#api-python-interop")[API / Python Interop].

The Nix workspace graph owns archive selection, feature validation and platform
filters. Explicit features outside the cached contexts are rejected before
building. Cargo validates its own feature selections. `just ci-update`
regenerates the workspace graph and scheduling configuration after context changes;
there is no separate development-test metadata file.

With Nix, the command builds or substitutes only that crate's archives and
runner dependencies, using the same feature contexts and runtime setup as CI.
Crates represented in multiple groups run each matching context, including the
integration/Python split with its prepared Python runtime. Test filters select
runtime execution without invalidating compilation artifacts. Tests execute on
every invocation, even when all compiled inputs and the runner are already cached.
An empty context is allowed when another context selects tests; an entirely
empty selection fails unless `--no-tests=pass` or `--no-tests=warn` is supplied.
`--profile ci_slow` also selects the slow lane's stack and worker settings.
Pass `-E` and `-P` separately from other short options; clusters such as `-vE`
are rejected so they cannot override mandatory platform exclusions.
The shell backend prints a fresh local directory under `target/test-reports/`,
with reports organized by group. Execution runs as the invoking user outside
the Nix builder, using a private writable source snapshot, home, cache and
temporary directory. `GAMMALOOP_WORKSPACE_ROOT`, `NEXTEST_WORKSPACE_ROOT` and
`INSTA_WORKSPACE_ROOT` point fixture and snapshot lookup at the relocated source
for that invocation. These are runner-owned overrides, not persistent shell
configuration; direct stub generators also honor the workspace-root overrides.
The source snapshot retains its CI group's context, including sibling
helper crates launched by documentation tests; compiled archives remain selected
per crate. A runtime resolves the license before replacing its home/config
directories: a nonempty `SYMBOLICA_LICENSE`, then `SYMBOLICA_LICENSE_SIGNED`,
then the caller's `${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license`.
The staged source is removed on exit, including generated failure state;
reports copied to the requested destination remain available.
`just test-nix CRATE` remains an explicit Nix-only command;
`just test-nix-all` selects the full CI aggregate.

The same runners are available as
`nix run .#test-CRATE -- REPORT_DIR NEXTEST_ARGS...`; for example:

```sh
nix run .#test-linnest -- target/test-reports/linnest -E 'test(layout)'
```

`just test-nix` uses the Nix shell backend. Runtime arguments are executable arguments,
not Nix evaluation inputs. `nix build .#test-CRATE` only prepares an executable;
it neither runs tests nor produces reports.
`just test-nix-runner` checks argument handling and archive reuse using Kurvst's
small test suite, plus documentation fixtures and the integration/Python split.
`just development-check` runs shell-command, wheel and publication regressions
without compiling workspace code.

`just build CRATE [CARGO_ARGS...]`, ordinary CLI builds, tests and wheels use the
shared profile. Explicit release recipes retain `release`; CLI variants retain
`ufo_support` and optional `python_abi`. `check`, `clippy`, `doc` and `run` use the same shell entry point.
`nix/development.nix` owns the native build/check/doc/wheel derivations; compatible
builds reuse CI library artifacts. Built executables are atomically copied into
`target/PROFILE/` (`target/debug/` for `dev`) with fresh modification times for
the CLI launcher's profile selection described in Common Commands. Each installed
binary retains its own GC root at `target/nix/bin/PROFILE/BINARY`; another build
cannot discard the closure of a binary it did not replace. The UFO-enabled CLI
wrapper pins Python 3.13, its matching UFO runtime and renderer data, while wheel
builds retain the selected packaging interpreter.

`just python-wheel gammaloop|community` builds a native wheel;
`just python-develop gammaloop|community` also installs it in an existing virtual
environment. `--python`, `VIRTUAL_ENV`, the matching `CONDA_PREFIX`, then root
`.venv` determine the installation interpreter; system installs fail. All three
backends install built wheels with uv, forcing reinstallation only of the
selected distribution; dependency resolution may still change other packages.
Python source changes require rebuilding and reinstalling;
deliberately editable workflows can use manual `maturin develop`.

The default `gammaloop` wheel builds the separate `gammaloop._gammaloop` extension, using the
installed compatible Symbolica Python package. Building a combined kernel is
optional and requires selecting `community` explicitly.

Community is that combined `symbolica` host under
`examples/notebooks/symbolica-host`, registering FeynKit, Linnet and Spenso as
`hepkit`, `graph` and `tensor`, using the shared profile. It does not include
GammaLoop, and these recipes do not offer a standalone `feynkit-py` wheel.
GammaLoop accepts `--abi3`. Existing
`build-api` variants retain their profile and ABI choices. See
#link("../../CONTRIBUTING.typ#native-python")[Native Python packaging] for
dependency provisioning and explicit `--no-deps` use with existing Community
environments; root lockfile synchronization can replace the combined host.
The default `development-wheel` flake package is a GammaLoop extension
derivation using the shared profile, even without an explicit wheel request.

Nix wheel outputs retain extracted native extensions so Nix records their
runtime library references. Installation keeps a GC root under
`ENV/.nix-python/DISTRIBUTION`; exported wheels keep
`OUTPUT/.nix-WHEEL_FILENAME`. Installation and root replacement are serialized;
failed installations retain their pending root and the previous environment root.
Before installation, the selected interpreter checks the Nix native payload's
loader compatibility. Full package imports must pass before the permanent root
is replaced. Builds on Darwin omit Linux wheel compatibility flags; runtime
compatibility on each deployment host remains necessary.
Retain these roots while using the wheels,
which depend on the native Nix store. Packaged initializers supply a default
`SYMBOLICA_DATA_DIRECTORY` for renderer data without overriding an explicit
value.

== Local CI preparation

On `itphlies`, run `just ci-checks-and-upload` before pushing CI-enabled work.
The configured Nix daemon downloads matching cached outputs as needed; the
command checks static inputs, builds missing runtime producers and uploads their
closures. Finish the upload before pushing so NixCI can reuse that
work. Keep substitution enabled rather than downloading the whole cache or
forcing rebuilds.

The local selection is `configuration.onlyBuild` in `nix/ci.nix` for the current
system: static checks, group Nextest producers, test dependency and binary
producers, runtime launchers and support/barrier outputs. Use the preview below
to inspect the selection for the current workspace graph. Publication follows
actual output references, including runtime closures and barriers;
licensed test results remain local execution output. Request-specific native
development builds and Python wheels are not default upload roots.

GitHub's shared producer job publishes a compact five-root closure:
`alphal00p-docs-cargo-artifacts`, `crane-ci-prebuild`, `gammaloop`,
`gammaloop-python-module` and the `gammaloop-nextest-binaries` check. Its 8 GiB
raw NAR limit is intentional. Self-contained full merged archives drop input
references, so this closure omits many intermediate archives selected by the
local upload command. A successful GitHub producer publication does not imply
that every local CI producer is available from the cache. Closure size and
coverage must be measured from realized outputs; a partial local store cannot
establish the full publication size or all omitted producers.

Use `just ci-upload-preview [scope]` to inspect these boundaries. The scope is
`local`, `github` or `both`, defaulting to `both`; for example,
`just ci-upload-preview github` inspects the compact publication. The report
lists comparable selected attributes, output paths and realized/missing status, and sums
deduplicated closure paths and raw NAR bytes for realized roots. Missing roots
make the measurement partial. Selected CI producers outside GitHub's realized
requisites are confirmed absent from its full closure only when all GitHub
roots are realized. Raw NAR bytes do not measure compressed transfer size.

The tracked-source preview has a different flake identity from raw `.#`
publication: revision, NAR hash and modification metadata can differ, and raw
directory filters can include untracked files. The documentation runner embeds
source metadata, so `publication.sourceRevisionDependent` lists it separately as
not comparable across source identities. It is excluded from comparable-root
realization counts and closure measurements. Reported bytes are a comparable-root
subtotal; the entire local upload size remains unknown when roots are excluded,
even when comparable coverage is complete. This does not change upload selection.

The preview evaluates offline with import-from-derivation disabled and queries
the local store. It never realizes or downloads selected artifact outputs,
uploads, or runs tests. The existing `develop.sh` entry may first bootstrap
pinned shell and Python inspection tools; that tooling bootstrap is the limited
exception to read-only inspection.

== Configuration and cache boundaries

#table(
  columns: (1fr, 1fr),
  table.header([File], [Responsibility]),
  [`flake.nix`], [Compose public packages, checks, apps and development shells],
  [`nix/rust-workspace.nix`], [Rust sources, features, build profiles and artifact reuse],
  [`nix/development.nix`], [Native build, check, Clippy, Rustdoc and wheel derivations],
  [`nix/documentation.nix`], [Documentation sources, reusable build artifacts and publication checks],
  [`nix/cargo-runtime.nix`], [Prepared Cargo environments for licensed doctest and documentation execution as the invoking user],
  [`nix/ci.nix`], [Test groups, required producers and NixCI scheduling],
  [`nix/ci-workspace-graph.json`, `nix-ci.nix`], [Generated workspace graph and self-contained NixCI configuration],
  [`scripts/develop.sh`], [Backend selection and tooling bootstrap],
  [`scripts/develop_nix.sh`], [Tracked source snapshots, Nix commands and native artifact retention],
  [`scripts/develop_dyndrv.sh`, `nix/dyndrv-bootstrap.nix`, `nix/dyndrv.nix`], [Pinned worker bootstrap, compiler environment and explicit dyndrv dispatch],
  [`scripts/python_wheel.py`], [Python interpreter selection, wheel packaging and installation],
  [`justfile`, `just/tests.just`], [Public native development commands],
  [`just/ci.just`], [Local checks, uploads, regeneration and reporting],
  [`.github/scripts/run-nix-ci-check.sh`, `.github/workflows/nix.yml`], [Realize check inputs, execute each runtime runner and export timing metadata],
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
the full documentation publication check, necessary producers and final success.
The documentation check validates generated references and examples, renders
latest and snapshot sites, checks their HTML, and exercises publication behavior.
It executes as the invoking user with a runtime license; reusable Cargo
artifacts are scheduled separately. `just ci-checks` selects the same runtime
suite; uploads select its runners and compiled producers. Public flake `checks`
exclude licensed runtime checks. `nix flake check` covers unlicensed build/static
checks, including packaging, standalone Rustdoc, the persistent Typst renderer
and WASM; it does not establish licensed test success.

`hydraJobs` and `hydraJobsMain` retain executable runner derivations at their
licensed attributes. NixCI invokes `nix-ci-check-*` packages through its run
jobs. The Actions matrix explicitly restores runtime rows omitted by Hestia's
cache filter and executes each runner after building or substituting it.
`.github/scripts/run-nix-ci-check.sh` separates build and runtime exit status;
either failure fails the job. A cached runner still runs its tests.

`just test-nix-all`, `just doctest-nix` and the documentation step of
`just ci-checks` execute `nix-ci-check-gammaloop-nextest`,
`nix-ci-check-gammaloop-doctest` and `nix-ci-check-alphal00p-docs`.
These wrappers print their local report directory and accept an explicit
destination as their first argument. The Nextest runner also accepts runtime
selection arguments after that destination.

`nix/cargo-runtime.nix` preserves Crane's vendoring, toolchain and artifact hooks
in a prepared, unlicensed environment. Its runner resolves the caller's license
and enters that environment with `nix develop` to unpack private sources and run
the licensed phase locally. It sets the same three workspace-root overrides
for that invocation, uses private home/cache/Cargo/temp directories, suppresses
environment dumps, and removes
temporary state on exit. License values are never evaluation inputs or
derivation attributes. This boundary does not remove older derivations or
outputs that already contained a license from existing stores or caches.

The generic documentation packages `alphal00p-docs-pages` and
`alphal00p-docs-snapshot-fixture` are executables accepting one absent or empty
output directory. For example, `nix run .#alphal00p-docs-pages -- target/docs-pages`
renders and validates pages locally. Compiled documentation inputs stay
cacheable; building the package prepares its runner rather than rendering pages.
The pages runner reads `ALPHAL00P_DOCS_CHANNEL` (`latest` by default) and requires
`ALPHAL00P_DOCS_SNAPSHOT_TAG` exactly when the channel is `snapshot`.

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

Test-only Rust modules belong in the `ownTest` source list in
`nix/rust-workspace.nix`. Production library inputs exclude those files; unit-test
producers include them. Native and WASM arithmetic, Python embedding and extension
linking, and the lean/full Typst renderer remain distinct configurations.

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

== Browser and widget coverage

`packages.x86_64-linux.nix-ci-check-gammaloop-nextest-browser` is a separate
licensed runtime executable, scheduled by the `gammaloop-nextest-browser` NixCI
launcher. Run it with
`nix run .#nix-ci-check-gammaloop-nextest-browser -- target/test-reports/browser`.
It reuses
per-crate nextest archives for Linnet rendering, Spynso and FeynKit. Nixpkgs
from `flake.lock` supplies Chromium, Python Playwright, AnyWidget, Marimo and typst-py;
there is no browser installation or package resolution at test time. The same
pinned STIX Two Math assets used by the documentation renderer enter through
`SYMBOLICA_DATA_DIRECTORY`. Browser font requests are fulfilled from that bundle.

The check runs the ignored real widget/Chromium regressions, the native MathML
producer/client round trip, tensor font/layout and library-grid checks, and
FeynKit CFF/LTD/scoped energy displays. Fixtures come from the archived native
owners, with disposable loopback servers and browser processes owned by the
runner. Fresh per-package JUnit identities select validators for the producers
that actually ran. The live MathML producer shares the caller's native test
selection while the runner coordinates its client. Browser runners reject
Nextest stress options before launching native tests.
Optional Typst-only tests also run in this dependency-complete group. Ordinary Spynso nextest runs include the paging harness; only
its live AnyWidget case skips when that optional package cannot import.
The embedded interpreter is the PyO3 compile interpreter, whose store closure
contains Python's standard library. The browser group adds its optional Python
packages through `PYTHONPATH`, without changing the compiled interpreter.

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

Each Nextest runtime copies available per-package JUnit reports with individual
case durations into its local report directory, including reports emitted by a
failed test run. GitHub Actions exports test names, statuses and
durations to a CSV artifact. Captured test output, expressions and raw logs
are excluded. A failed preparation or runtime still fails the job, even when
there is no report or timing CSV to upload. Reports are execution output outside
the Nix store.
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
not the commit. Runners retain the archives and runtime inputs needed to execute
their selected tests. The local upload command selects static checks, runtime
runners and compiled producers; runtime results are not cacheable Nix checks.
Building or substituting these outputs is preparation only. Every NixCI run
executes the selected runtime suite using those inputs and requires its success.
Test groups are scheduled and reported independently.

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
comparisons: the existing physical workload and nine scalar route cases, plus
two cases that build and benchmark every evaluator method and one thermal
case. Scalar route cases
share their point corpus and settings. The runner renders `scalar.toml` with
each workload's graph and route settings; the other cases use fixed cards.
The exact rendered card and declared inputs, including DOT graphs and shared
scalar model and restriction files, are hashed before measurement. Both
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
Graph enumeration is outside this timing because every workload imports a fixed
DOT input. GL15 retains the original top-pentagon topology, projector, loop basis
and numerical corpus independently of generator numbering. Its import and
integrand construction are still measured on both revisions.
The companion `benchmarks/performance/gg-hhh-topology.toml` card supports a
separate local topology-generation comparison. It enumerates the complete
ungrouped gg → hhh process with the same explicit projector. Compare fresh CLI
processes in five alternating baseline/candidate pairs, and verify the exported
particle/vertex topology multiset with fixed external labels and overall factors
outside the timed process. This diagnostic is separate from the thirteen-case
CI gate and does not certify symbolic numerator equivalence by topology alone.
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
Every physical result must be finite, nonzero and agree between revisions.
Scalar results use relative tolerance 10⁻⁹ with no absolute allowance, because
an absolute floor suitable for GL15 could mask changes to their smaller values.
Each topology/route comparison, including any confirmation run, has a
four-minute budget. The collection can take longer than four minutes.
Building the two CLI binaries is outside that budget. Exceeding it yields
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
SHAs. The separate `Performance main` workflow measures every main push, including documentation-only merges, against the previous main SHA. These captures run independently of PR and merge-group gates. Ordinary feature-branch pushes
do not trigger this workflow. Manual benchmark dispatch remains available.

A relevant new commit while `final-review` remains set must rerun the gate.
Returning to development means removing that label; PR-number concurrency
cancels superseded work. Compiled CLI packages are reused from Nix caches, but
the benchmark always measures fresh states. All thirteen cases run sequentially on
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
source filters, run `just ci-update`, and commit `nix/ci-workspace-graph.json`
and `nix-ci.nix`. There is no separate development test metadata file. FeynKit
needs its eighth group, nine extracted crates, embedded fixtures, optional Python
dependencies, boundary checks and Python stub checks. Its extracted libraries
must remain independent of GammaLoop and its private workspace-hack dependency.
The shared FeynKit branch was not changed by the experiments.

The `ci-cache-base` input pins a compatible compiler-state seed. Keep that commit
reachable in published history. Refresh deliberately after a green revision with
`just ci-cache-base FULL_COMMIT_SHA`; rebasing or cleaning history alone is not a
reason to invalidate it. See the #link("nix-crane-cache-reuse.typ")[cache-reuse audit]
for earlier investigation and implementation details.

The optional `ci-compiler-state` package collects compiler state for a compatible
baseline; neither ordinary local uploads nor the compact GitHub producer job
select it. Publish that baseline once after its revision is publishable, then
pin it with `just ci-cache-base FULL_COMMIT_SHA`. Do not publish every commit's
compiler state as part of the ordinary upload cycle.

== Measured results

These measurements describe earlier implementations and revisions, not validation
of the current FeynKit Nix port. The historical warm implementation could reuse
cached test-result checks. Current runners always execute licensed tests; their
runtime cost and correctness require fresh evidence for the candidate revision.

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
- Historical test execution is classified as `executed`, `reused` or `unknown`.
  Fetching a runner wrapper proves neither execution nor test-result reuse;
  current licensed runners must provide fresh execution evidence. Nextest
  `failed` includes `timedOut`;
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

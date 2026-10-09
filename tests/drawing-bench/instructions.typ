= Drawing benchmark P0

Python 3.10+ and Typst 0.15.x are required. This harness adapts the supplied
`run.py`: fresh processes, one untimed warmup per variant, case-by-case
interleaving and at least three rounds. It does not implement a fuel/prototype
harness; none exists in the supplied bundle.

```sh
python3 tests/drawing-bench/run.py --dry-run --out target/drawing-bench/setup
python3 tests/drawing-bench/run.py --rounds 3 --png 150
python3 tests/drawing-bench/run.py --cases higgs-sunset-particles --out target/drawing-bench/single
python3 tests/drawing-bench/run.py --variants current current-stock --stock-package-path /path/to/stock-packages --out target/drawing-bench/paired
```

Output directories must be new and under `target/drawing-bench`. No references
are updated. The deterministic archive is SHA256-verified before safe extraction;
all 24 fixtures remain byte-for-byte unchanged. A dry run assembles the same
current sources and metadata without compiling. Inspect extracted cases and
licenses in the output's `fixtures/` directory.

== Current versus historical

Normal runs copy CURRENT `crates/linnest/typst/src`,
`crates/kurvst/typst/src`, their checked-in WASM and
`assets/embedded/drawing/templates`. The harness stages the unmodified CeTZ
0.5.2 archive bundled for offline Python rendering, plus MiTeX and oxifmt,
in the run's private package directory. No editable CeTZ fork is maintained
by this benchmark. Ensure these WASM assets match your source build before
measuring; their hashes are recorded rather than implicitly rebuilding them.

`current-stock` changes only the package path. Supply upstream CeTZ 0.5.2,
MiTeX 0.2.6 and oxifmt 1.0.0 there. It is an independent stock distribution
comparison; both variants use the same current renderer. Never substitute the
old scalar overlay to make this mode succeed: its old optimize arguments are
incompatible with consolidated core, and it would measure the wrong sources.
Neither variant uses an overlay or a second rendering pipeline.

The original external bundle may be measured separately with its own `run.py`,
and named HISTORICAL/FROZEN snapshot, never CURRENT. This harness neither
invokes nor depends on it.

== Captured reconciled baseline

`references.tar` contains the 24 current-source PNGs at 150 ppi, all 72 raw
timing samples, case medians, source/package/environment metadata, and the
paired comparison with the supplied snapshot. SHA256:
`633b5b736b31d9f1378cf9b4108c3785a036e0d874558b90b469297b3be3142c`.

The sum of three-round case medians is 15.010 s on the recorded M4 Pro host
with Typst 0.15.1. All 24 PNGs are byte-identical to that host's supplied
snapshot renders. This certifies no raster drift at 150 ppi in the reconciled
baseline, not cross-platform SVG identity or future arrowhead implementation.

```sh
mkdir -p target/drawing-bench/reference-images
tar -xf tests/drawing-bench/references.tar -C target/drawing-bench/reference-images
```

== Measurement and review

`timings.csv` stores every measured sample in execution order, in milliseconds;
`medians.csv` reports each case's median, sum of medians and last/first ratio.
The first requested variant is the paired baseline. Keep diagram pairs
(momenta and particles/notebook) together when reporting grouped regressions.
PNG export is optional and untimed. Wall time includes process startup,
compilation, SVG export, file output and process exit; assembly and hash
inventory are outside the timing boundary. Warmups exclude the initial cold
disk reads, not subsequent process startup.

`metadata.json` records revision, actual assembled source/package/fixture
identities, executable hash/version, host/OS, sample count, order, fonts and
output hashes. Working-tree changes are represented by source hashes, not
assumed to match HEAD. System fonts are ignored and `TYPST_FONT_PATHS` removed.

Hashes identify artifacts, not approval. Historical Linux SVG references are
retained unchanged, not treated as pixel-identity acceptance tests.
User-approved policy allows bounded arrow/label differences plus explicit
image review. Review paired images at identical resolution, recording case,
reviewer, displacement/shape bounds, accepted differences and rejected
overlaps/clipping in a review record alongside the run. No numerical tolerance
has been approved here: do not invent one. Link targets and hit-area coverage
must also be inspected in SVG/browser output; PNG equality alone cannot certify
them. The harness always marks visual review pending.

== Provenance

Imported read-only from `/Users/lcnbr/Downloads/drawing-bench`.
Supplied README: feynkit `ec21ce9a`, 2026-10-06; drawing sources and vendored
packages unchanged since `78c33d71`. Cases captured 2026-10-01 from the former
CeTZ route, subsequently adjusted only to import MiTeX directly and use
root-relative imports.

Historical EPYC 9754 Linux run on 2026-10-08: vendored 36.1 s, stock 95.4 s,
2.64×. Earlier lighter-load measurements: 26.7 s / 70.0 s. These are
host-specific frozen references, not new Mac measurements or current-source
claims. The separately captured reconciled PNGs are stored in `references.tar`.

Linnest/Kurvst licenses and Linnest's ec-layout/clarabel notices are preserved
inside the archive under their original `tree/crates/...` paths. Existing
repository packages retain their own CeTZ LGPL-3.0, MiTeX and oxifmt licenses.

== Focused native Python acceptance

`native.py` is independent of `run.py`, the Typst corpus, package selection,
and P7 packaging. Use an installed community Symbolica host providing
`symbolica.community.hepkit`, and the installed `linnet` extension, on Python
3.10+. Keep any user license in the process environment, never in artifacts.
Pause other benchmarks and CPU-heavy builds before measuring.

```sh
python3.13 tests/drawing-bench/native.py
# A repeat requires a new output directory:
python3.13 tests/drawing-bench/native.py --out target/drawing-bench/native-acceptance-repeat
```

The fixed scalar cube-minus-edge DOT avoids diagram generation. The checked-in
`crates/feynkit-py/tests/fixtures/scalars_2p_3p.json` model must validate its
four loops, eleven internal edges, eight cubic interaction vertices and two
external edges. Render configuration is `layouts: {steps: 100}`, momentum
arrows on, particle and momentum labels off, no title. This has empty label
pages and uses the native `Scene::render` branch, not Typst typesetting.
FeynKit currently accepts a JSON configuration dictionary, not a
`linnet.RenderConfig` instance; the dictionary uses the equivalent layout
and template options without changing either API.

Ten untimed warmups precede twenty-one `perf_counter_ns` samples. Only the
installed Python `diagram.render(config=...)` call is timed; imports, model
loading, validation, SVG parsing, hashing and writes are excluded. The reported
median is compared with an absolute upper target of 100 ms. This is not a
nonregression certificate: no genuine old same-host Python baseline is
available. Rust unit tests and labelled/Typst renders are not timing substitutes.

Each new directory under `target/drawing-bench/native-acceptance` (or the
explicit `--out`) retains input DOT/configuration, arrows-on/off SVGs, all
sample SVGs and `metadata.json`: raw nanosecond samples, median, host/toolchain,
loaded module paths/hashes, fixture/source hashes, revision and dirty-diff hash.
Installed binary hashes identify the actual extension; the harness does not
assume it was built from current dirty sources. Preserve its build log alongside
the run to establish that relationship.

Untimed XML checks require a finite positive viewbox, thirteen nondegenerate
painted momentum-head paths (the default open stroked V, not a filled triangle),
no corresponding heads with arrows off, and stable node/edge/half-edge identities.
Interactive links require parseable details and positive hit rectangles; any
remaining fragment links must resolve. These checks are not visual approval of
overlaps, clipping, shaft contacts or hit-area coverage; review the SVGs separately.
Exit status is 0 for a validated median <= 100 ms, 1 for an exceeded target,
and 2 for a blocked import/render/assertion. A blocked run records no invented
median. Exception text and environment values are deliberately not serialized
because license errors may contain private data.

=== Local host setup and execution evidence

The existing canonical community host is
`examples/notebooks/symbolica-host/pyproject.toml`; its Maturin manifest is
`examples/notebooks/symbolica-host/Cargo.toml`. Build/install only into a
workspace-local environment, for example:

```sh
python3.13 -m venv target/drawing-bench/native-host-build/venv
CARGO_TARGET_DIR="$PWD/target/drawing-bench/native-host-build/cargo" \
  maturin build --locked --release \
  --manifest-path examples/notebooks/symbolica-host/Cargo.toml \
  --interpreter "$PWD/target/drawing-bench/native-host-build/venv/bin/python" \
  --out target/drawing-bench/native-host-build/wheels
# Install the resulting host wheel and a matching linnet wheel into that venv,
# then run native.py using its bin/python. Do not install globally.
```

Initially this macOS host's Python 3.13 had no installed Symbolica community host.
The bounded offline locked build stopped because the host lockfile needed
resolution after path-dependency changes. An authorized offline unlocked
attempt then stopped before compilation: public dependency `jiff-static
v0.2.35` is not cached. Offline resolution also downgraded existing dependencies,
so the original lockfile was restored rather than retaining unrelated churn.
Logs, exit statuses, original/resolved lockfiles and the resolution diff are
under `target/drawing-bench/native-host-build`. No Python render performance
acceptance or nonregression result is established by these setup attempts.

A subsequent normal online build used a generated copy of the canonical host
under `target/drawing-bench/native-host-build/host`, with path dependencies
rebound to this workspace. Public downloads succeeded, including `jiff-static`.
Its generated lockfile added 95 packages (and updated `rust-embed-utils` to
match its newly resolved family); no checked-in lockfile was changed.
The first release build reached the native dependency stack before its 150 s
bound; a single resumed locked build reached FeynKit generation/CFF/tensor
compilation before its 180 s bound. Both exited 124 without producing a host
wheel. These timeouts were incomplete setup attempts, not evidence of an
unavailable dependency.
`build-online.log`, `build-online-resume.log` and their exit-status files retain
the process evidence. The generated manifest, lock and partial Cargo cache can
be reused with a suitably longer bounded build; do not treat them as an installed
host or as benchmark results.

The initial `python3.13 tests/drawing-bench/native.py` invocation on this host
exited 2 at the first import (`linnet` was also not installed). Its
`target/drawing-bench/native-acceptance/metadata.json` records `status: blocked`,
without timing samples or a median. Ruff formatting/lint and Python syntax
validation passed; isolated SVG-inspector checks exercised real open-V geometry,
degenerate marks, empty hitboxes and unresolved fragments, not Python rendering.

A final sufficiently long locked release invocation completed the real
community host in 55.14 s, reusing the partial cache. A matching standalone
Linnet Maturin release build completed in 4 min 14 s. Both wheels were installed
with `--no-deps --no-index` exclusively into
`target/drawing-bench/native-host-build/venv`. Their canonical sources were not
modified; the host uses community Symbolica revision `942bd2c0` and the pinned
Symbolica Typst plugin revision from its canonical manifest. Full build logs,
exit statuses and source/manifest/lock identities are retained as
`build-final.log`, `build-linnet.log` and `build-identity.json`.

```sh
target/drawing-bench/native-host-build/venv/bin/python tests/drawing-bench/native.py \
  --out target/drawing-bench/native-acceptance-installed
```

That real installed invocation exited 2 before rendering/timing. Symbolica
reported that the user license key format is outdated and must be renewed at
`https://symbolica.io/license/`. No key or environment values are recorded.
This is now the actual blocker: provide a valid renewed user license through
the process environment and rerun with a fresh `--out` directory. Installed
module paths/hashes are present in
`target/drawing-bench/native-acceptance-installed/metadata.json`; its
`status: blocked` and `error_type: RuntimeError` have no samples or median.
No absolute 100 ms acceptance or Python nonregression claim is established.
The build succeeds; fixed-DOT validation and native rendering still require
the valid license. Parent Typst corpus timing remains independent.

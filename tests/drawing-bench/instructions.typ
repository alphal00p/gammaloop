= Drawing benchmark P0

Python 3.11+ and Typst 0.15.x are required. This harness adapts the supplied
`run.py`: fresh processes, one untimed warmup, case-by-case order and at least
three rounds. It does not implement a fuel/prototype harness; none exists in
the supplied bundle.

```sh
python3 tests/drawing-bench/run.py --dry-run --out target/drawing-bench/setup
python3 tests/drawing-bench/run.py --rounds 3 --png 150
python3 tests/drawing-bench/run.py --cases higgs-sunset-particles --out target/drawing-bench/single
python3 tests/drawing-bench/run.py --package-path /path/to/typst/packages --out target/drawing-bench/stock
python3 tests/drawing-bench/run.py --inspection-targets --out target/drawing-bench/targets
python3 tests/drawing-bench/run.py --label-report --out target/drawing-bench/labels
python3 tests/drawing-bench/run.py --label-report \
  --label-baseline target/drawing-bench/labels/label-report.json --out target/drawing-bench/labels-new
```

Output directories must be new and under `target/drawing-bench`. No references
are updated. The deterministic archive is SHA256-verified before safe extraction;
all 24 fixtures remain byte-for-byte unchanged. A dry run assembles the same
current sources and metadata without compiling. Inspect extracted cases in the
output's `fixtures/` directory.

== Current sources and stock packages

Runs copy the CURRENT `crates/linnest/typst/src`, `crates/kurvst/typst/src`,
their checked-in WASM and `assets/embedded/drawing/templates`. Ensure these
WASM assets match your source build before measuring; their hashes are
recorded rather than implicitly rebuilding them.

Packages come from `--package-path`, a Typst package directory with stock
CeTZ 0.5.2, MiTeX 0.2.6 and oxifmt 1.0.0 from Typst Universe, the only packages
the sources import. It defaults to `TYPST_PACKAGE_CACHE_PATH`, which the
development shell points at its pinned package set. The harness checks the
three packages before compiling and removes `TYPST_PACKAGE_CACHE_PATH` from the
compile environment. No CeTZ fork or overlay is part of this benchmark. The original external bundle may be measured separately with its
own `run.py`, and named HISTORICAL/FROZEN snapshot, never CURRENT.

== Captured reconciled baseline

`references.tar` contains the 24 current-source PNGs at 150 ppi, all 72 raw
timing samples, case medians, source/package/environment metadata, and the
paired comparison with the supplied snapshot. SHA256:
`633b5b736b31d9f1378cf9b4108c3785a036e0d874558b90b469297b3be3142c`. These
references are the P0 renders with the formerly vendored, patched CeTZ 0.5.1;
stock CeTZ renders are compared against them by image review.

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
`medians.csv` reports each case's median and their sum. Compare revisions by
running the harness from each checkout on the same host, interleaved, and keep
diagram pairs (momenta and particles/notebook) together when reporting grouped
regressions. Runs are static by default, but revisions whose Linnest `draw` has
no `inspection-targets` parameter, and the recorded reference medians, compiled
every case with its link targets; time this checkout with `--inspection-targets`
when comparing against them (`metadata.json` records the flag as
`inspection_targets`). PNG export is optional and untimed. Wall time includes process
startup, compilation, SVG export, file output and process exit; assembly and
hash inventory are outside the timing boundary. Warmups exclude the initial
cold disk reads, not subsequent process startup.

`metadata.json` records the checked-out commit (from jj in a jj workspace,
otherwise from Git; `null` for a tree outside any checkout), actual assembled
source/package/fixture identities,
executable hash/version, host/OS, sample count, order, fonts and output hashes.
Working-tree changes are represented by source hashes, not assumed to match
HEAD. System fonts are ignored and `TYPST_FONT_PATHS` removed.

Hashes identify artifacts, not approval. Historical Linux SVG references are
retained unchanged, not treated as pixel-identity acceptance tests.
User-approved policy allows bounded arrow/label differences plus explicit
image review. Review paired images at identical resolution, recording case,
reviewer, displacement/shape bounds, accepted differences and rejected
overlaps/clipping in a review record alongside the run. No numerical tolerance
has been approved here: do not invent one. The cases are static documents, so
their renders carry no link targets. `--inspection-targets` compiles them with
`--input linnest-inspection-targets=true`, as the interactive viewers do; inspect
link targets and hit-area coverage in that SVG/browser output, since PNG equality
alone cannot certify them. The targets are invisible, so both modes give the same
PNGs. The harness always marks visual review pending.

`--label-report` adds an automated, untimed label-overlap check after timing:
it evaluates each case once more with `--input linnest-label-report=true` and
queries Linnest's `<linnest-label-report>` metadata. The shared label placement
computes it with the geometry its search scores: painted edge lines, stroke
radius included, clipped exactly against each chosen label's text box; and
other labels' text, nodes, self-loop bends and momentum arrows whose boxes
share a positive area with the label's text or its own arrow. Arrow and node
boxes are rectangles around their shapes, so a listed box overlap need not be
visible. `label-report.json` keeps every drawing's records; `label-report.csv`
lists per case the labels, those whose text meets a painted edge, overlapping
label pairs, each crossing as `e<label edge>>e<crossed edge>:<length inside the
box>`, each obstacle as `e<label edge>>n<node>` or `e<label edge>>e<edge>-loop`,
and each momentum arrow as `e<label edge>>e<arrow edge>`. With
`--label-baseline`, the run fails when a case gains a (label edge, crossed
edge, obstacle or arrow edge) contact or an overlapping label pair relative to
an earlier `label-report.json`. The report describes the chosen candidates only;
whether a clear candidate existed needs Linnest's `debug-label-collisions`
overlay or the placement's candidate costs.

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

Reproduced on 2026-10-09 on itphlies (AMD EPYC 9754, Linux, Typst 0.15.1)
from the consolidated Kurvst sources and the then vendored packages: all 24
PNGs at 150 ppi are byte-identical to `references.tar`. The three-round sum of
case medians was 27.2 s at a shared-host load average near 58; it is not
comparable with the M4 Pro figure and is recorded only as context for paired
runs on this host.

The supplied archive also held copies of the Linnest/Kurvst licenses,
Linnest's clarabel notice and the notice of the since removed ec-layout plugin
under `tree/`. The archive contains none of that code, so they are not retained
(see the manifest). Stock packages carry their own CeTZ LGPL-3.0, MiTeX and
oxifmt licenses.

== Focused native Python acceptance

`native.py` is independent of `run.py`, the Typst corpus and package selection.
Install the Symbolica Community host wheel (`examples/notebooks/symbolica-host`),
which provides `symbolica.community.graph` and `symbolica.community.hepkit`, into
a dedicated environment on Python 3.10+, and run the script with that
environment's interpreter. Keep any user license in the process environment,
never in artifacts. Pause other benchmarks and CPU-heavy builds before
measuring.

```sh
python tests/drawing-bench/native.py
# A repeat requires a new output directory:
python tests/drawing-bench/native.py --out target/drawing-bench/native-acceptance-repeat
```

The fixed scalar cube-minus-edge DOT avoids diagram generation. The checked-in
`crates/feynkit-py/tests/fixtures/scalars_2p_3p.json` model must validate its
four loops, eleven internal edges, eight cubic interaction vertices and two
external edges. Rendering uses `RenderSettings(layouts=LayoutSettings(impred_steps=100))`
and `DiagramStyle` with momentum arrows on, particle and momentum labels off
and no title. This has empty label pages and uses the native `Scene::render`
branch, not Typst typesetting.

Ten untimed warmups precede twenty-one `perf_counter_ns` samples. Only the
installed Python `diagram.render(config=..., style=...).to_svg(interactive=True)` call is timed;
imports, model loading, validation, SVG parsing, hashing and writes are
excluded. The reported median is compared with an absolute upper target of
100 ms. This is not a nonregression certificate; compare installed wheels of
both revisions on the same host for that.

Each new directory under `target/drawing-bench/native-acceptance` (or the
explicit `--out`) retains input DOT/configuration, arrows-on/off SVGs, all
sample SVGs and `metadata.json`: raw nanosecond samples, median,
host/toolchain, loaded module paths/hashes, fixture/source hashes and revision.
Installed binary hashes identify the actual extension; the harness does not
assume it was built from current dirty sources. Preserve its build log
alongside the run to establish that relationship.

Untimed XML checks require a finite positive viewbox, thirteen nondegenerate
painted momentum-head paths (the default open stroked V, not a filled
triangle), no corresponding heads with arrows off, and stable
node/edge/half-edge identities. Interactive links require parseable details
and positive hit rectangles; any remaining fragment links must resolve. These
checks are not visual approval of overlaps, clipping, shaft contacts or
hit-area coverage; review the SVGs separately. Exit status is 0 for a validated
median <= 100 ms, 1 for an exceeded target, and 2 for a blocked
import/render/assertion. A blocked run records no invented median. Exception
text and environment values are deliberately not serialized because license
errors may contain private data.

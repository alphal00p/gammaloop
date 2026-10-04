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
`assets/embedded/drawing/templates`. They reuse the repository package path
`crates/linnet-py/vendor/typst-packages`; no renderer tree or CeTZ is duplicated
in this benchmark. Ensure these WASM assets match your source build before
measuring; their hashes are recorded rather than implicitly rebuilding them.

`current-stock` changes only the package path. Supply upstream CeTZ 0.5.1,
MiTeX 0.2.6 and oxifmt 1.0.0 there. Until P7 removes vendored-only calls,
compilation may fail explicitly with a Typst diagnostic. Never substitute the
old scalar overlay to make this mode succeed: its old optimize arguments are
incompatible with consolidated core, and it would measure the wrong sources.
After migration this same mode needs no overlay or second rendering pipeline.

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

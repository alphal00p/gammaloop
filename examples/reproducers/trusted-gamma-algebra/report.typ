= Trusted gamma algebra and alias domains

The 2026-09-29 comparison measures Idenso's public `simplify_gamma` against
commit `0bfc7003`. Both builds use Symbolica
`6a96c9d77b217ea394770524c119e0613d4d9ef5` and the same frozen dependency
artifacts. Idenso and Spenso use optimization level 2 with debug assertions;
this is a native Rust comparison, not Python or FORM timing.

The change retains validation facts on the shared symbolic tensor, including
when exposing alias roots and definitions. Trusted gamma identities publish
their established interface directly. New homogeneous trace sums read one
representative term with the existing fast structure reader. The trace emitter
also carries the input's intrinsic-normalization proof, avoiding another walk
of every generated sum. User normalizers and unresolved positional interfaces
retain checked boundaries.

== Measurements

Times are milliseconds per warm public call, including result destruction.
The table excludes initial input admission and final alias resolution; both
implementations return the same aliased result abstraction. Admission-inclusive
variants are saved separately. Initial evaluation warms lazy trace kernels.
Three process blocks each contain five samples; the reported value is the
median of the three block medians. Each process is pinned to CPU 12 on a shared
host, with interleaved baseline/candidate blocks.

#table(
  columns: (2.5fr, 1fr, 1fr, 1fr),
  [Case], [Before (ms)], [After (ms)], [Speedup],
  [4D trace, 8 free indices], [24.025], [3.429], [7.01×],
  [4D trace, 12 free indices], [443.408], [60.285], [7.36×],
  [4D axial trace, length 12], [1278.923], [793.741], [1.61×],
  [D-dimensional trace, 8 free indices], [20.574], [2.867], [7.18×],
  [D-dimensional trace, 10 free indices], [77.122], [9.635], [8.00×],
  [D-dimensional trace, alternating p/q, length 12], [1.616], [0.453], [3.57×],
  [Open chain, repeated indices, length 6], [1.128], [0.816], [1.38×],
  [1 alias domain], [0.781], [0.310], [2.52×],
  [16 alias domains], [8.745], [3.677], [2.38×],
  [64 alias domains], [39.555], [18.655], [2.12×],
)

All 17 case/mode variants preserve the exact resolved expression and logical
interface, and are idempotent under a second public simplification. Those
checks run outside the timer. Alias cases contain the indicated number of
reachable scalar gamma definitions plus one disconnected definition. The
ordinary/D-dimensional, axial and open-chain cases use the same default gamma
settings. D-dimensional cases retain a spin trace unit of four.

The repository regressions separately assert zero checked-inference and
scope-validation calls for admitted gamma traces and certified alias-domain
passes. They retain independent checked-interface oracles, malformed-input
rejections, callback-induced rank-loss cases, logical AUTO order and typed
zeros. The benchmark is a same-input regression comparison, not a new
independent Clifford-algebra proof.

Final qualification passed 809 Idenso unit tests (22 skipped), all 78 Spenso
unit tests, strict Clippy for both crates' libraries and tests, and formatting.
The existing exhaustive gamma5-position test exceeded the standard four-minute
timeout; a local per-test override allowed twenty minutes without changing
assertions. The complete Idenso test phase then passed in 645 seconds.
Both changed Typst documents compile, and the portable linker passes Ruff.

Axial length 12 remains the slowest measured case. This harness does not
profile its residual algebra or alias bookkeeping, so the remaining time is
not assigned to an unmeasured internal phase. Full production pipelines,
component execution, cold startup and FORM comparisons are outside this
measurement.

== Reproduction

#link("timings.json")[`timings.json`] records all samples, block medians,
build features, artifact hashes and the protocol. The measurement driver is
#link("driver.rs")[`driver.rs`]. Its checked-in formatting differs from the
measured source only by rustfmt, recorded in
#link("source-formatting.json")[`source-formatting.json`]. The final production
source differs from the compiled benchmark snapshot only by a line wrap;
the other recorded source differences are test-only.

To avoid rebuilding Symbolica or changing feature unification,
#link("link-driver.py")[`link-driver.py`] links against explicit compatible
Idenso, Spenso, Symbolica and serde-json libraries. Build both implementations
with the same Rust toolchain and dependency family, then supply their exact
library paths and every required dependency directory:

```sh
python link-driver.py --idenso /path/libidenso.rlib \
  --spenso /path/libspenso.rlib --symbolica /path/libsymbolica.rlib \
  --serde-json /path/libserde_json.rlib --dependency-dir /path/deps \
  --output /tmp/gamma-baseline
```

Repeat for the candidate, using its Idenso/Spenso libraries and the same
Symbolica dependencies. The linker records its command and input hashes.
Run qualification before timing:

```sh
/tmp/gamma-baseline /tmp/gamma-reference check
/tmp/gamma-candidate /tmp/gamma-check check /tmp/gamma-reference
taskset -c 12 /tmp/gamma-baseline /tmp/gamma-before measure
taskset -c 12 /tmp/gamma-candidate /tmp/gamma-after measure /tmp/gamma-reference
```

Use separate output directories for each block and retain stdout JSON lines.
The saved block order and timing/admission distinction are in `timings.json`.

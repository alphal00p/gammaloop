# GL638 build10 backend timing

Measured on 14 September 2026 using one authenticated executable and the same
512 complete Samples, each repeated three times: 1536 draws per worker setting.
This is the previous **1000 GeV, seven-channel** benchmark (six cuts plus Cut1
joint H/Z), mu_r=m_uv=91.188 GeV, all 936 orientations and six cuts. It is not a
600 GeV production estimate or a recommendation to omit explicit soft coverage.
State, physical settings, sampling maps and inputs were unchanged between arms.

The native source/host checks are documented in
[GL638_NATIVE_PRECISION_ACCEPTANCE.md](GL638_NATIVE_PRECISION_ACCEPTANCE.md).
Build10 permits the requested O3 level through the shared process evaluator
activation route; the backend-only standalone loading/export routes still use
fixed O2. Symbolica/SymJIT dependency pins were unchanged for this comparison.

## Measured components

S is parameterization time, E is time inside evaluators across all physical
calls, and P−E is the remaining integrand evaluation time. Entries are mean
**elapsed worker milliseconds per draw**, not exclusive CPU time. All 13,824
timing evaluations remained in Double precision; rescue counts cannot explain
the concurrency differences.

| Backend | Workers | S ms/draw | E ms/draw | P−E ms/draw | Global wall s | Global draws/s |
|---|---:|---:|---:|---:|---:|---:|
| Eager | 1 | 7.90 | 32.85 | 3.44 | 68.99 | 22.3 |
| Eager | 20 | 8.67 | 76.04 | 9.21 | 15.59 | 98.5 |
| Eager | 50 | 9.95 | 102.79 | 19.23 | 35.62 | 43.1 |
| SymJIT O2 compressed | 1 | 7.91 | 20.20 | 3.81 | 50.11 | 30.7 |
| SymJIT O2 compressed | 20 | 8.43 | 25.56 | 6.06 | 12.77 | 120.2 |
| SymJIT O2 compressed | 50 | 9.90 | 32.18 | 11.61 | 29.31 | 52.4 |
| SymJIT O3 compressed | 1 | 7.98 | 19.45 | 3.80 | 49.04 | 31.3 |
| SymJIT O3 compressed | 20 | 8.11 | 21.45 | 6.67 | 11.99 | 128.1 |
| SymJIT O3 compressed | 50 | 10.00 | 28.73 | 15.62 | 32.97 | 46.6 |

Global wall includes clone destruction; it excludes clone construction/warmup,
state loading and JIT activation. It therefore does not provide a production
throughput estimate for a long iteration. Summed worker durations and parallel
global wall are different quantities and must not be added together.

O3 evaluator time is about 16% below compressed O2 at 20 workers and 11% below it
at 50 workers in this deck. At 50 workers, however, S+P is 53.69 ms/draw for O2 and
54.34 ms/draw for O3: the lower evaluator time is offset by other measured time.
This single short run supplies no uncertainty estimate or statistical ranking
of that small whole-call difference. Compressed O3 is a reasonable candidate
for the physics screen; O2 remains a useful control. Both substantially reduce
measured evaluator time relative to Eager.

## Why 50-worker global throughput appears worse

| Backend | Workers | Longest recorded worker chain s | Global wall minus that chain s | Excluded clone/warmup s |
|---|---:|---:|---:|---:|
| Eager | 20 | 8.58 | 7.02 | 33.81 |
| Eager | 50 | 5.69 | 29.93 | 78.03 |
| SymJIT O2 compressed | 20 | 3.62 | 9.15 | 33.80 |
| SymJIT O2 compressed | 50 | 2.53 | 26.78 | 76.59 |
| SymJIT O3 compressed | 20 | 3.27 | 8.72 | 34.38 |
| SymJIT O3 compressed | 50 | 2.66 | 30.30 | 78.78 |

The scratch driver consumes each worker clone inside its parallel closure and
returns only batch results. Destruction of that clone therefore occurs after
its last batch timer stops but before the global timer stops. With Eager, the
longest evaluation chain is **shorter** at 50 workers (5.69 s) than at 20 (8.58 s),
yet the global wall is longer because the uninstrumented remainder grows from
7.02 s to 29.93 s. O2 and O3 show the same large remainder.

Clone teardown is a source-supported explanation for this discrepancy, rather
than evidence that 50 physical evaluations per unit time must be slower than 20.
The remainder also includes task scheduling and result collection; without
teardown timestamps it cannot all be assigned quantitatively to deallocation.
The data do show longer elapsed evaluator calls at higher concurrency, but do
not isolate L3-cache saturation from memory bandwidth, scheduling, warmup or
concurrent shared-host activity.

Production owns clones differently: `integrate/mod.rs` constructs
`CoreIterationState` values **inside each integration iteration**, retains them
across chunks, merges their runtime state, then drops them at iteration end.
Thus clone cost is amortized per iteration, not only once per hour. Calibrate
production from actual iteration wall times and avoid extrapolating these
short-deck global rates as a constant cost per draw. No implementation changes
or further optimization were made for this audit.

## Repetition and numerical checks

Evaluator milliseconds per draw by repetition:

| Backend | Workers | Repetitions 0 / 1 / 2 |
|---|---:|---|
| Eager | 20 | 94.91 / 76.17 / 57.04 |
| Eager | 50 | 136.26 / 101.11 / 70.98 |
| SymJIT O2 compressed | 20 | 26.39 / 25.09 / 25.19 |
| SymJIT O2 compressed | 50 | 36.29 / 32.35 / 27.88 |
| SymJIT O3 compressed | 20 | 21.52 / 21.53 / 21.31 |
| SymJIT O3 compressed | 50 | 33.09 / 28.86 / 24.23 |

Repetitions are serial within each worker, without a global barrier between
repetitions. Two initial points warm each clone, but Eager in particular is not
stationary across the measured repetitions. They are not three independent
throughput trials. One-worker Eager is stable at 32.83–32.87 ms inside evaluators.

The existing three-arm artifact analyzer passed 44,712 checks and 13,312 paired
numerical comparisons, with zero failed configured-budget comparisons. Each
arm's prior native-control audit also passed. The hard complete real controls
agree with forced-Arb values; individual effectively-zero Cut2 real residues
retain pure-relative misses (absolute differences no larger than 1.27e-19).

The timed deck uses its historical complex-norm stability settings. Its 18
own-component misses all repeat **source 64** across O2/O3, worker counts and
repetitions: real estimates approximately −1.18327e-39 versus Eager −4.48363e-40,
an absolute difference 7.34910e-40. These misses remain explicitly recorded; the benchmark does not certify every
real component separately. The 600 GeV physics controls use componentwise
stability and provide separate evidence.

A subsequent 600 GeV, nine-channel componentwise control run also passed all
11 retained Samples with O3. Complete real and absolute-real estimators, and
precision/accuracy histories, matched O2 exactly in serialized output (nine
Double and two Quad evaluations). This panel used no absolute comparison floor
and did not integrate. Its independent audit is
`/tmp/gl638-final-physics-screen/component-controls-o2-o3-build10-independent-audit.json`,
SHA256 `0f762ddd0d7c450917dca77b9614e072ecd9bf0607c051f2d8761ca23e026e89`.

## Evidence

- Analyzer: `/tmp/gl638-symjit-next/timing-analysis-build10-three-arm/analysis.json`.
- Arm records: `/tmp/gl638-symjit-next/results-{eager,symjit_o2_compress_on,symjit_o3_compress_on}-build10/summary.json`.
- Executed requests: `request-timing-eager-build10.json` and `request-timing-jit-build10.json` in the same directory; they differ only in status/provenance hashes.
- Driver: `gate_timing_build10.rs`, source SHA256
  `33d3caf26da6e89270b8ff93612b067488d1c9f8749a2269367b5dcf42a17cbd`;
  executable SHA256 `8b0c8a74b96f2aa4a762769410543bc9c1d18abbe3f7287cbbb41480d782e7f3`.
- JIT activation: 11.14 s O2, 10.48 s O3. Complete arm elapsed times were 315.38 s
  and 318.44 s, including loading, controls and setup. The wrapper reports peaks
  around 146.7 GiB; its O3 peak is cumulative across the two child runs, so it
  cannot establish a precise per-backend RAM difference.

These results provide compiler and execution-cost evidence only. They contain
no new cross-section central value, MC error, variance improvement or global
weight-bound claim.

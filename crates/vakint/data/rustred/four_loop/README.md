# Four-loop RustRed candidate programs

These assets supply the four-loop branch of the opt-in
`EvaluationMethod::RustRed` / `EvaluationOrder::rustred_only()` backend.
They do not change Vakint's default evaluation order or its FORM-backed modes.

Each parent (H, FG, BMW and X) has:

- a `.csv` ordered physical/auxiliary momentum descriptor;
- a `.toml` topology-generic RustRed family-generation input;
- a `.candidates.toml.gz` saved parametric candidate program;
- a `.rrcat` exact map of its declared finite terminals onto FMFT's PR basis.

Labels identify data files, not engine dispatch. Vakint selects a program from
the defining parent and simultaneous routing already retained by its topology
matcher. It does not rematch the graph. The descriptor, loaded family and
terminal catalog must agree exactly.

## Runtime contract

Each program is decompressed and loaded once, on first use, through RustRed's
candidate loader. A shared native RustRed applier then handles scalar-numerator
lowering, guard selection, strictly descending rule application and memoization.
No rule-generation campaign or FORM executable is invoked. Vakint consumes the
result and reuses its pure-Rust/Symbolica FMFT master finalizer. FeynKit's normal
tensor prepass permits a fully FORM-less evaluation chain.

These are **candidate programs, not `ClosedArtifact`s**. The loader does not
assert independent regenerated-source replay or unlimited family closure.
Runtime guards, descent checks and the explicit terminal set remain enforced;
missing rules or unknown terminal values fail rather than falling back to FORM
or inventing masters. Numerical comparison tests establish the tested inputs,
not completeness for every numerator rank. On 2026-09-19, the public backend
passed all fifteen original four-loop numerical references and all sixteen
expanded-propagator/pinch comparisons against FMFT, using an invalid FORM path
for FeynKit and RustRed. Full-process times were 132.59 s and 208.69 s
respectively, including both peers and program loading, not scalar-only
timings. See the [acceptance commands and evidence](../../../tests/experimental_rustred_4l/README.md).

`RustRedEvaluationOptions { substitute_masters: false }` leaves the raw PR master
basis unexpanded. With substitution enabled, the existing finite FMFT expansion
tables are used. The `.rrcat` records are exact PR expressions, **not new
20,000-digit master evaluations**. Requesting more precision than a surviving
master/constant source contains emits the existing warning through Vakint's
logger; it does not increase source accuracy. Missing Laurent orders are errors.

## Offline reproduction

The GammaLoop runtime pin `13bc9474` remains compatible with these assets.
For the offline example commands below, use RustRed `91751eb2` or a compatible
later revision: that helper amendment raises its collection-entry allowance
for the large four-loop bundles without changing the runtime loader or schema.
From that RustRed checkout, with the Symbolica license supplied in the
environment:

```sh
cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  generate /path/to/vakint/data/rustred/four_loop/h.toml 9 6 default \
  /path/to/new/h.candidates.toml /path/to/new/h.report.toml

cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  verify 10 /path/to/new/h.candidates.toml \
  /path/to/vakint/data/rustred/four_loop/h.rrcat

gzip -n -9 -c /path/to/new/h.candidates.toml > /path/to/new/h.candidates.toml.gz
```

The helper requires new output paths. The zero-based nonpositive index list is
`9` for H/X and `8,9` for FG/BMW; these are auxiliary scalar-product coordinates,
not physical propagators. The worker count above is an explicit example and can
be changed. Generation uses ordinary RustRed IBPs and no FMFT-derived rules.
The fresh-process `verify` operation checks the family/catalog binding and
declared terminal keys; it is not a symbolic closure certificate.

FMFT is needed only when preparing or independently revalidating the offline PR
catalog and in the separate oracle test lane. It is not a dependency of runtime
RustRed evaluation. Catalogs must be refreshed after any correction to the
oracle's numerator-routing conventions; candidate programs need not be
regenerated for a catalog-only correction. The 2026-09-19 audit re-evaluated
all 1,155 entries in 199.29 s, correcting four FG projections after fixing its
signed numerator routing. H, BMW and X had no differences; a subsequent FG
audit checked all 145 entries with zero differences in 23.62 s. These are exact
PR-projection checks, not new high-precision numerical master evaluations.

## Public scalar timing probe

The ignored public test measures four dotted parents, four expanded propagator
numerators and a factorized four-tadpole. It prepares each input once through
the normal matcher and FeynKit, then supplies the identical scalar expression
to `EvaluationMethod::FMFT` and `EvaluationMethod::RustRed`. Numerical agreement
is required for every observation; no rule generation or experimental reducer
injection occurs.

Run this exact test in a fresh process, with `SYMBOLICA_LICENSE` supplied in the
environment and a workspace-local `TMPDIR`:

```sh
VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH=/absolute/path/to/form \
RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
CARGO_PROFILE_RELEASE_LTO=off CARGO_BUILD_JOBS=4 \
  cargo test --release --locked -p vakint --test rustred_four_loop_tests \
  public_timing::representative_public_scalar_timings -- \
  --exact --ignored --nocapture --test-threads=1
```

`PUBLIC_SCALAR_TIMING` records contain case, phase, backend, repetition, wall
nanoseconds, optional CPU seconds and parent RSS snapshots before/after the
call. CPU includes waited-for children; snapshots are not phase-specific peak
RSS and omit a concurrently live FORM child's memory. Measure whole-process
peak RSS separately. Snapshots accumulate already loaded parent programs;
CPU resolution is the platform clock tick reported in the metadata.
Compilation is outside the timing records.

First-parent-use observations include lazy program loading and first point
application; they do not isolate loader cost. Five alternating paired repeated
calls retain process-local caches. The public API does not expose point-cache
clearing or rule counters, so these repeats must not be described as cold-point
timings. An initial call for a subsequent input may also reuse subproblems
visited by earlier inputs. Tensor preparation and final numerical comparison
are outside the scalar timer; normal public settings validation, FORM startup and I/O, scalar
reduction and master substitution remain inside it. Any matching or dispatch
performed internally by `evaluate_integral` is also included: this is a public
backend measurement, not a bare `CandidateReducer` kernel timing.

### Measured public-backend snapshot, 19 September 2026

All nine inputs passed the initial comparison and five paired repetitions at
relative tolerance `1e-20`, with zero uncertainty allowance. The release build
used optimization level 3 with LTO disabled, on an AMD EPYC 9754 host, pinned to
physical cores 88–93 with one scalar caller and nested thread pools capped at
one. Shared-host contention was not otherwise eliminated.

| Input | First RustRed call (s) | Repeated RustRed median (ms) | Repeated FMFT median (ms) | FMFT / RustRed |
| --- | ---: | ---: | ---: | ---: |
| H, first propagator cubed | 35.353 | 106.624 | 1483.568 | 13.91 |
| H, expanded D7 numerator | 0.021 | 21.055 | 138.442 | 6.58 |
| FG, first propagator cubed | 11.450 | 86.678 | 225.339 | 2.60 |
| FG, expanded D7 numerator | 1.661 | 32.252 | 147.771 | 4.58 |
| BMW, first propagator cubed | 20.907 | 79.288 | 716.347 | 9.03 |
| BMW, expanded D7 numerator | 0.055 | 41.740 | 140.989 | 3.38 |
| X, first propagator cubed | 127.140 | 156.553 | 4717.944 | 30.14 |
| X, expanded D7 numerator | 2.683 | 78.877 | 169.894 | 2.15 |
| Factorized four-tadpole | 0.016 | 14.816 | 134.338 | 9.07 |

Each cubed-parent first call includes that parent's lazy load and uncached
application. Subsequent inputs can reuse earlier subproblems. In particular,
the 127.140 s X observation is **not a loader-only timing**. Repeated calls are
faster than FMFT on this finite matrix, but first-use costs remain substantial.
Even with the parent program already loaded, initial FG/X expanded-D7 calls
took 1.661/2.683 s, versus FMFT's 0.149/0.194 s. The repeated-call advantage
therefore must not be generalized to previously unseen targets.

The complete benchmark process took 251.79 s wall, 231.27 s user CPU and 18.63 s
system CPU; peak RSS was 14,686,292 KiB, with no swaps. This process includes
all nine inputs, both backends, initialization and untimed comparisons. CPU
observations have 10 ms granularity. Per-call RSS records are snapshots, not
independent backend peaks. Raw observations and the derived summary are kept
in the local RustRed workspace's ignored
`TMP/four-loop-public-benchmark.{log,time}` and
`TMP/four-loop-public-benchmark-summary.tsv`; the test above reproduces the
workload without those files.

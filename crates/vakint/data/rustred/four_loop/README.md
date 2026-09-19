# Four-loop RustRed candidate programs

These assets supply the four-loop branch of the opt-in
`EvaluationMethod::RustRed` / `EvaluationOrder::rustred_only()` backend.
They do not change Vakint's default evaluation order or its FORM-backed modes.

Each parent (H, FG, BMW and X) has:

- a `.csv` ordered physical/auxiliary momentum descriptor;
- a `.toml` topology-generic RustRed family-generation input;
- a `.candidates.rrbin.gz` saved native-binary parametric candidate program;
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

The binary programs contain Symbolica's native rational-polynomial atoms and
shared state, native family geometry and structural rule/source records. They
are trusted, build-time embedded data for the pinned RustRed/Symbolica stack;
the native decoder is not an untrusted-file parser. Decompression and framing
have explicit size limits. Native transport does not change rule applicability,
the terminal basis or certification status.

These are **candidate programs, not `ClosedArtifact`s**. The loader does not
assert independent regenerated-source replay or unlimited family closure.
Runtime guards, descent checks and the explicit terminal set remain enforced;
missing rules or unknown terminal values fail rather than falling back to FORM
or inventing masters. Numerical comparison tests establish the tested inputs,
not completeness for every numerator rank. On 2026-09-19, the public backend
passed all fifteen original four-loop numerical references and all sixteen
expanded-propagator/pinch comparisons against FMFT after the native-binary
migration, using an invalid FORM path for FeynKit and RustRed. Full-process
times were 32.45 s and 104.98 s respectively, including both peers and program
loading, not scalar-only timings. These suites ran concurrently on disjoint
CPU sets. Three fixture checks and three loader/catalog/raw-master unit checks
also passed. See the [acceptance commands and evidence](../../../tests/experimental_rustred_4l/README.md).

The later atomic K1/K3/K6 certified-V6 migration pins RustRed `d6718733` and
leaves these four candidate files and all master catalogs unchanged. The full
fifteen-reference and sixteen-pinch public gates pass again on that pin, at
34.42 s and 104.55 s full-process wall time, respectively. An independent
auditor ran them sequentially on CPUs 80–85 with nested pools capped at one;
the through-three-loop selection ran concurrently on disjoint CPUs. These are
correctness-run timings, not a new paired performance comparison.

Migration preserved all 59,636 rules and 1,155 declared terminal keys. Every
saved coefficient, source, guard and family component was compared exactly
offline; no rules were regenerated. The four compressed programs total
20,511,348 bytes, versus 29,947,045 bytes for the former text programs, using
the same deterministic `gzip -n -9` command.

`RustRedEvaluationOptions { substitute_masters: false }` leaves the raw PR master
basis unexpanded. With substitution enabled, the existing finite FMFT expansion
tables are used. The `.rrcat` records are exact PR expressions, **not new
20,000-digit master evaluations**. Requesting more precision than a surviving
master/constant source contains emits the existing warning through Vakint's
logger; it does not increase source accuracy. Missing Laurent orders are errors.

## Offline reproduction

Use the matching RustRed runtime pin `d6718733bee4d20f4554e9d694b288ddaae60475`
for the offline example commands below, with the Symbolica license supplied
in the environment. These commands generate fresh programs when intentionally
requested; the shipped migration itself converted the saved programs without
generating new rules.

```sh
cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  generate /path/to/vakint/data/rustred/four_loop/h.toml 9 6 default \
  /path/to/new/h.candidates.rrbin /path/to/new/h.report.toml

cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  verify 10 /path/to/new/h.candidates.rrbin \
  /path/to/vakint/data/rustred/four_loop/h.rrcat

gzip -n -9 -c /path/to/new/h.candidates.rrbin > /path/to/new/h.candidates.rrbin.gz
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

### Native-binary public-backend snapshot, 19 September 2026

All nine inputs passed the initial comparison and five paired repetitions at
relative tolerance `1e-20`, with zero uncertainty allowance. The release build
used optimization level 3 with LTO disabled, on an AMD EPYC 9754 host, pinned to
physical cores 88–93 with one scalar caller and nested thread pools capped at
one. Shared-host contention was not otherwise eliminated.

| Input | First RustRed call (s) | Repeated RustRed median (ms) | Repeated FMFT median (ms) | FMFT / RustRed |
| --- | ---: | ---: | ---: | ---: |
| H, first propagator cubed | 15.962 | 105.526 | 1669.229 | 15.82 |
| H, expanded D7 numerator | 0.023 | 21.514 | 143.815 | 6.68 |
| FG, first propagator cubed | 1.643 | 83.668 | 244.900 | 2.93 |
| FG, expanded D7 numerator | 1.610 | 32.444 | 152.148 | 4.69 |
| BMW, first propagator cubed | 6.288 | 75.208 | 716.279 | 9.52 |
| BMW, expanded D7 numerator | 0.055 | 41.284 | 136.986 | 3.32 |
| X, first propagator cubed | 82.864 | 150.665 | 4662.828 | 30.95 |
| X, expanded D7 numerator | 2.934 | 84.236 | 177.525 | 2.11 |
| Factorized four-tadpole | 0.015 | 15.136 | 132.988 | 8.79 |

Each cubed-parent first call includes that parent's lazy load and uncached
application. Subsequent inputs can reuse earlier subproblems. In particular,
the 82.864 s X observation is **not a loader-only timing**. Repeated calls are
faster than FMFT on this finite matrix, but first-use costs remain substantial.
Even with the parent program already loaded, initial FG/X expanded-D7 calls
took 1.610/2.934 s, versus FMFT's 0.169/0.240 s. The repeated-call advantage
therefore must not be generalized to previously unseen targets.

The complete benchmark process took 167.33 s wall, 161.29 s user CPU and 4.85 s
system CPU; peak RSS was 2,863,816 KiB, with no swaps. This process includes
all nine inputs, both backends, initialization and untimed comparisons. CPU
observations have 10 ms granularity. Per-call RSS records are snapshots, not
independent backend peaks. Raw observations and the derived summary are kept
in the local RustRed workspace's ignored
`TMP/gamma-native-migration-20260919/public-timing.{log,time}`; the test above
reproduces the workload without those files.

The previous same-workload text-program snapshot took 251.79 s overall and
14,686,292 KiB peak RSS; its first H/FG/BMW/X calls took
35.353/11.450/20.907/127.140 s. These separate shared-host runs are diagnostic
comparisons, not a controlled statistical speedup claim. Native transport
removes text-loading overhead; it does not change the application algorithm or
make uncached reductions as fast as the repeated-cache timings.

# Four-loop RustRed candidate programs

These assets supply the four-loop branch of the opt-in
`EvaluationMethod::RustRed` / `EvaluationOrder::rustred_only()` backend.
They do not change Vakint's default evaluation order or its FORM-backed modes.

Each parent (H, FG, BMW and X) has:

- a `.csv` ordered physical/auxiliary momentum descriptor;
- a `.toml` topology-generic RustRed family-generation input;
- a `.candidates.rrbin.gz` saved native-binary parametric candidate program;
- a `.rrcat.bin` native exact map of its declared finite terminals onto FMFT's PR basis.

Labels identify data files, not engine dispatch. Vakint selects a program from
the defining parent and simultaneous routing already retained by its topology
matcher. It does not rematch the graph. The descriptor, loaded family and
terminal catalog must agree exactly.

## Runtime contract

Each program is decompressed and loaded once, on first use, through RustRed's
candidate loader. At that boundary, RustRed also prepares and verifies a
same-family vacuum terminal-equivalence plan once. A shared native RustRed applier then handles scalar-numerator
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

The terminal-value files use RustRed's distinct `TerminalValues` envelope and
deduplicate exact Symbolica Atoms with a shared native state. Vakint delegates
their encode/decode to RustRed; runtime loading performs no expression-text
parsing. Their family fingerprint, arity, declared coverage and all 1,155 raw
integral keys are unchanged. These catalogs carry no IBP or closure authority.
The consumer still checks equality with the candidate's declared terminal set.

The catalog-only migration compared every value exactly against the saved text,
then repeated that comparison in a fresh process with unrelated Symbolica state:

| Catalog | Keys | Distinct exact values | Former text bytes | Native bytes |
| --- | ---: | ---: | ---: | ---: |
| H | 386 | 23 | 25,632 | 8,845 |
| FG | 145 | 18 | 10,309 | 5,654 |
| BMW | 179 | 20 | 13,514 | 6,391 |
| X | 445 | 21 | 33,565 | 9,417 |
| Total | 1,155 | 26 across parents | 83,020 | 30,307 |

This 63.5% reduction is a transport/dictionary improvement, not a smaller
physical master basis. No FMFT calculations or IBP generation were repeated.
Native atom bytes may depend on ambient Symbolica IDs; imported exact values
and their namespaces do not. Earlier text catalogs are deliberately unsupported.

The catalog migration pins RustRed `d51721b6`. Its default-feature public backend
again passes all fifteen numerical references and all sixteen expanded-numerator/
pinch comparisons, with forbidden FORM paths in both native stages. The two
whole-process correctness gates took 33.04 s and 115.23 s respectively, including
the separate FMFT oracle, with peak RSS of 1,990,636 and 2,689,952 KiB. These are
shared-host validation observations, not a new paired performance benchmark.
The unchanged 83-test through-three-loop selection, seven catalog/loader unit
checks and three four-loop fixture checks also pass. Terminal aliases were
disabled for that catalog-only milestone; the later activation is described below.

## Exact vacuum terminal normalization

The `2b50267c` runtime pin explicitly steers
`TerminalAliasPlan::vacuum_parametric_equivalences` with its default preparation
bounds. RustRed verifies equality of the full restricted first Symanzik
polynomial U, including its scale, under a power-preserving parameter
permutation. For supported positive-power unit-mass vacuum keys this proves a
unit-coefficient integral equality; it does not assert an integer loop-momentum
map. Native Symbolica graph canonicalization proposes permutations and exact
polynomial comparison verifies them. Equivalent terminals coalesce inside the
existing memoized applier, before ancestor coefficient accumulation. Vakint
implements neither the proof nor new coefficient arithmetic.

| Parent | Raw catalog keys (unchanged) | Prior routing representatives | U representatives | Verified U aliases |
| --- | ---: | ---: | ---: | ---: |
| H | 386 | 176 | 52 | 334 |
| FG | 145 | 53 | 26 | 119 |
| BMW | 179 | 68 | 37 | 142 |
| X | 445 | 208 | 64 | 381 |
| Total | 1,155 | 505 | 179 | 976 |

These are family-local representatives, **not 179 independent masters**. They
comprise 74 positive-power keys and all 105 unchanged negative-index keys.
Unsupported shapes and exhausted preparation bounds retain the original keys;
none of these four saved parents reached a bound. All 976 aliases were checked
against the exact saved catalog expressions in the core gate. The raw-key
coverage check, saved files, guard checks, descent, source conditions and
unresolved-leaf errors remain unchanged. No programs or master values are
regenerated, and the plan adds no family-closure certificate.

The production loader opts into `from_reducer_with_terminal_aliases` on a fresh
owner. The existing experimental `from_reducer` constructor retains its original
behavior, including accepting a populated cache. Installing aliases requires
an explicitly empty cache and never silently clears caller data. Vakint's
evaluation options, default method order and FORM-backed methods are unchanged.

Plan preparation adds first-use work, including for inputs already equal to a
terminal. It is amortized across later reductions of that parent; the public
benchmark below includes it in the first-parent call rather than presenting
isolated core application times as whole-backend timings.

The independently executed U-activation four-loop gates pass all fifteen
original numerical references and all sixteen expanded-numerator/pinch pairs.
The unchanged 83-test selection through three loops, nine focused constructor/
catalog/loader checks and three four-loop fixture checks also pass.
The two four-loop processes took 22.56/46.90 s wall time, with peak RSS
1,859,076/2,092,664 KiB; these correctness runs include the separate FMFT oracle.
FeynKit and RustRed use forbidden FORM paths. No numerical tolerance or expected
value changed. The previous routing-only activation (`f91c47ab`) passed the same
gates at 24.71/47.88 s; those validation observations are not a controlled timing
comparison. The separately matched public benchmark appears below.

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
tables are used. The `.rrcat.bin` records are exact PR expressions, **not new
20,000-digit master evaluations**. Requesting more precision than a surviving
master/constant source contains emits the existing warning through Vakint's
logger; it does not increase source accuracy. Missing Laurent orders are errors.

## Offline reproduction

Use the matching RustRed runtime pin `2b50267c5b1df10e4cd2798d7863d3f8c1b47383`
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
  /path/to/vakint/data/rustred/four_loop/h.rrcat.bin

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

### Matched terminal-alias comparison, 19 September 2026

Both the frozen no-alias runtime (`d51721b6`) and the alias-enabled runtime
(`f91c47ab`) passed all nine inputs, each with an initial comparison and five
paired repetitions: 108 timed calls and 54 numerical comparisons per process.
Relative tolerance is `1e-20` with zero uncertainty allowance. Both release
builds use optimization level 3 with LTO disabled, on an AMD EPYC 9754 host,
affinity 88–93, one scalar caller and nested pools capped at one. Shared-host
contention was not eliminated; other validation/build work used disjoint CPU
sets. These separate process runs are diagnostics, not confidence bounds.

Times shown as before → after compare the same saved programs and public test
inputs. The FMFT column gives five-call medians from both runs, so movement in
the reference backend is visible too.

| Input | First RustRed call (s), before → after | Warm RustRed median (ms), before → after | Warm FMFT median (ms), before → after |
| --- | ---: | ---: | ---: |
| H, first propagator cubed | 14.859 → 7.155 | 103.061 → 90.010 | 1536.791 → 1490.813 |
| H, expanded D7 numerator | 0.021 → 0.021 | 20.953 → 20.826 | 139.285 → 134.571 |
| FG, first propagator cubed | 1.244 → 1.339 | 83.896 → 80.593 | 234.677 → 220.815 |
| FG, expanded D7 numerator | 1.715 → 0.782 | 32.275 → 30.931 | 147.154 → 145.719 |
| BMW, first propagator cubed | 7.326 → 3.623 | 75.680 → 69.522 | 721.270 → 715.237 |
| BMW, expanded D7 numerator | 0.054 → 0.032 | 41.319 → 27.482 | 139.454 → 135.494 |
| X, first propagator cubed | 84.167 → 30.990 | 150.270 → 126.704 | 4710.187 → 4726.422 |
| X, expanded D7 numerator | 2.707 → 0.872 | 76.313 → 54.745 | 181.094 → 163.144 |
| Factorized four-tadpole | 0.015 → 0.015 | 14.616 → 15.312 | 130.457 → 142.288 |

The observed H/X first-call improvements are 2.08×/2.72×, including lazy load
and new plan preparation. They are **not loader-only improvements**. FG's
first cubed-line call is slightly slower, and the terminal-heavy factorized
control gains nothing: preparation and scheduling still cost time. The public
dotted probes have power three; separate core-only diagnostics using power two
are a different workload and must not be substituted for this table.

Repeated alias-enabled calls are faster than FMFT on this finite matrix, but
first use is still slower than FMFT for every cubed parent. For example, initial
H/X FMFT calls took 1.502/4.734 s versus RustRed's 7.155/30.990 s. Initial FG/X
expanded-D7 calls also remain slower at 0.782/0.872 s versus 0.144/0.190 s for
FMFT. The warm-cache advantage is not a bound for unseen targets.

| Whole benchmark process | No aliases | Aliases enabled |
| --- | ---: | ---: |
| Wall time | 164.83 s | 96.35 s |
| User + system CPU | 158.66 + 4.83 s | 92.19 + 3.53 s |
| Peak RSS | 2,864,780 KiB | 2,218,696 KiB |

These process totals include both backends, initialization and untimed
comparisons. RSS includes loaded programs, caches and temporary values, not
only cached coefficient payload. Per-call CPU has 10 ms granularity; RSS
snapshots are not independent backend peaks. Neither run regenerated rules or
master catalogs. The workspace-only Cargo cache was initialized between the
runs; copying dependencies, compilation and download time are excluded.

Raw observations, exact executable hashes, phase tables and an identical
before/after driver are retained in the RustRed workspace's ignored
`TMP/gamma-terminal-alias-rollout.4Zi7ej/` as `baseline-d517.*`, `aliases-f91.*`
and `run-benchmark.sh`. The public test above reproduces the workload without
those local files. The earlier native-I/O-only snapshot remains in Git history
and `TMP/gamma-native-migration-20260919/`; it is not the baseline used here.

### Matched U-plan comparison, 19 September 2026

This later comparison uses the **routing-enabled** `f91c47ab` runtime as its
control, not the earlier no-alias `d51721b6` runtime. The new pin is `2b50267c`.
A fresh process was run for each, using the identical nine inputs, unchanged
saved programs/catalogs and the same five alternating paired repetitions.
Both pass all 54 numerical comparisons (108 timed calls). Affinity, release
profile, single scalar caller, nested pool caps and timing boundaries are
identical to the preceding experiment. Compiler and independent validation
processes occupied disjoint CPU sets on the same shared host; the table is
diagnostic evidence, not a statistical speedup guarantee.

| Input | First RustRed call (s), routing → U | Warm RustRed median (ms), routing → U | Warm FMFT median (ms), routing → U |
| --- | ---: | ---: | ---: |
| H, first propagator cubed | 7.177 → 6.263 | 90.211 → 82.422 | 1502.284 → 1476.941 |
| H, expanded D7 numerator | 0.022 → 0.021 | 20.902 → 20.370 | 148.245 → 132.621 |
| FG, first propagator cubed | 1.429 → 1.042 | 80.646 → 76.979 | 229.253 → 221.330 |
| FG, expanded D7 numerator | 0.833 → 0.765 | 31.417 → 31.417 | 173.022 → 144.075 |
| BMW, first propagator cubed | 3.967 → 3.247 | 73.561 → 63.989 | 772.049 → 713.649 |
| BMW, expanded D7 numerator | 0.032 → 0.030 | 29.689 → 26.132 | 167.386 → 136.908 |
| X, first propagator cubed | 31.272 → 29.450 | 125.005 → 116.638 | 4678.632 → 4798.671 |
| X, expanded D7 numerator | 0.857 → 0.634 | 54.329 → 51.662 | 161.097 → 162.069 |
| Factorized four-tadpole | 0.015 → 0.015 | 15.591 → 14.810 | 133.482 → 129.176 |

The first H/X calls improve by about 1.15×/1.06× in these observations, including
lazy loading and plan preparation. FG's repeated expanded-D7 cost is effectively
unchanged. Repeated RustRed calls are faster than FMFT on all nine tested
inputs, but initial cubed-parent calls still cost more than FMFT: H/X
6.263/29.450 s versus 1.485/5.052 s. These warm observations do not predict
unseen-point costs. The smaller representative set does not imply a proportional
speedup because loading, nonterminal recursion and coefficient arithmetic remain.

| Whole benchmark process | Routing plan | U plan |
| --- | ---: | ---: |
| Wall time | 98.78 s | 95.83 s |
| User + system CPU | 94.32 + 3.81 s | 91.82 + 3.40 s |
| Peak RSS | 2,215,440 KiB | 2,161,216 KiB |

These process totals include **both backends**, initialization and untimed
comparisons; they are not isolated RustRed application times. Neither run
regenerated rules or master values. Raw records, executable hashes, exact
commands and independently recomputed median tables are retained in
`TMP/gamma-terminal-u-rollout.oFoQwr/` as `baseline-f91.*` and
`u-2b50267.*`. The unchanged public command above reproduces the workload.

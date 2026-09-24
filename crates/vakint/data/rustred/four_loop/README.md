# Four-loop RustRed candidate programs

These assets supply the four-loop branch of the opt-in
`EvaluationMethod::RustRed` / `EvaluationOrder::rustred_only()` backend.
They do not change Vakint's default evaluation order or its FORM-backed modes.

The current runtime pin is `7b22f5bda441588f2e437b15fa0872461c483708`.
The September 24 equation-to-package refresh supplies 59,509 unrestricted
candidate rules in 20,381,288 compressed bytes. All 1,155 raw terminals,
74 catalog outputs and normalization sidecars are retained. See the
[current producer recipe and measured validation](../../../rustred_package_refresh.typ).
The migration and timing sections below retain their historical measurements;
their 59,636-rule payload is not the current generated package.

Each parent (H, FG, BMW and X) has:

- a `.csv` ordered physical/auxiliary momentum descriptor;
- a `.toml` topology-generic RustRed family-generation input;
- a `.candidates.rrbin.gz` saved native-binary parametric candidate program;
- a `.rrcat.bin.gz` compressed native exact map of its normalized output terminals onto FMFT's PR basis;
- a `.rrnorm.bin` native weighted terminal-normalization sidecar.

Labels identify data files, not engine dispatch. Vakint selects a program from
the defining parent and simultaneous routing already retained by its topology
matcher. It does not rematch the graph. The descriptor, loaded family and
terminal catalog must agree exactly.

## Runtime contract

Each program and its small value catalog are decompressed and loaded once,
on first use. RustRed owns candidate loading and the native value codec.
At that boundary, RustRed also decodes and independently
rebuilds the same-family weighted vacuum terminal-normalization plan once.
A shared native RustRed applier then handles scalar-numerator
lowering, guard selection, strictly descending rule application and memoization.
No rule-generation campaign or FORM executable is invoked. Vakint consumes the
result and reuses its pure-Rust/Symbolica FMFT master finalizer. FeynKit's normal
tensor prepass permits a fully FORM-less evaluation chain.

The binary programs contain Symbolica's native rational-polynomial atoms and
shared state, native family geometry and structural rule/source records. They
are trusted, build-time embedded data for the pinned RustRed/Symbolica stack;
the native decoder is not an untrusted-file parser. Decompression and framing
have explicit size limits. Candidate-program transport alone does not change
rule applicability, declared raw terminals or certification status; the separate
normalization sidecar defines the effective output convention.

The terminal-value files use RustRed's distinct `TerminalValues` envelope and
deduplicate exact Symbolica Atoms with a shared native state. Vakint delegates
their encode/decode to RustRed; runtime loading performs no expression-text
parsing. Their family fingerprint, arity, declared coverage and retained exact
values are unchanged. Only the 74 possible normalized outputs require catalog
entries; all 1,155 raw keys and their relations remain in the separate program
and normalization sidecar. These catalogs carry no IBP or closure authority.
The consumer checks exact equality with RustRed's installed output-terminal
set, rejecting missing values and unused raw entries. No legacy catalog fallback
exists. Single-member gzip transport rejects truncated/corrupt/trailing data
and applies a 64-KiB input/decompressed catalog bound before native import.

## Compressed output-only catalogs

The value-only migration neither generates rules nor changes the terminal
normalizer. Every retained value is copied as an exact native Atom. All 1,155
previous raw projections agree exactly with the 74-value catalog composed with
the unchanged sidecars (1,260 coefficients checked). Original master precision
and numerical tolerances are unaffected.

| Parent | Raw declarations retained in program | Catalog outputs | Previous native bytes | Output-only native bytes | Shipped gzip bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| H | 386 | 22 | 8,845 | 5,777 | 929 |
| FG | 145 | 16 | 5,654 | 5,115 | 821 |
| BMW | 179 | 17 | 6,391 | 5,431 | 861 |
| X | 445 | 19 | 9,417 | 5,609 | 897 |
| Total | 1,155 | 74 | 30,307 | 21,932 | 3,508 |

The shipped value catalogs are 88.4% smaller than the previous uncompressed
native files. This is a catalog-storage improvement, not an 88.4% reduction of
the much larger candidate programs or a new master-basis reduction. Compression
has zero timestamp and no filename metadata. The runtime retains the existing
per-program lazy owner and RustRed memoized rule applier.

With the same compressor applied to the original raw catalogs, they would total
6,157 bytes; pruning therefore reduces even that compressed baseline by 43.0%.
Without compression the native subset alone is 27.6% smaller. The four large
compressed candidate programs are unchanged at 20,511,348 bytes, as are the
30,795 bytes of normalization sidecars. The combined binary payload changes
from 20,572,450 to 20,545,651 bytes; pruning values does not shrink rule storage.

Three fresh processes per family/transport give the following median times from
file read through native import, excluding process launch (milliseconds):

| Parent | Previous native catalog | Output-only native catalog | Shipped compressed catalog |
| --- | ---: | ---: | ---: |
| H | 0.742 | 0.587 | 0.608 |
| FG | 0.535 | 0.545 | 0.550 |
| BMW | 0.540 | 0.543 | 0.562 |
| X | 0.795 | 0.564 | 0.592 |

These are small shared-host observations, not a statistical speed guarantee or
cold-filesystem test. Imports start in fresh Symbolica processes, and gzip adds
approximately 0.067 ms of transport work. Other compiles used disjoint CPUs.
Some smaller inputs get slower; catalog compression is primarily a storage
improvement and is not claimed to accelerate the full integral reduction.

The 2026-09-20 release gate passes the unchanged 83-case through-three-loop
selection, fifteen focused catalog/normalization/transport checks, three
four-loop fixture checks, all fifteen original numerical references and all
sixteen expanded-numerator/pinch pairs. The public scalar harness additionally
passes 54 paired numerical comparisons (108 measured calls across nine inputs).
The RustRed/FeynKit branches retain invalid FORM paths; only the separate FMFT
oracle receives an executable. Input precision, nonunit scales, tolerances,
default evaluation order and existing backends are unchanged. The optional
offline-catalog audit target also compiles with `experimental-rustred` enabled.
Focused filters overlap the inventory and are not extra physical inputs.

The reference/pinch processes took 60.74/82.94 seconds wall time and peaked at
1,859,864/2,115,584 KiB RSS. These shared-host correctness runs used CPU 94 with
concurrent compiler contention, not an isolated performance comparison. The
focused unit process moved from CPU 94 to 92 near completion. The timing
harness also reloads raw rule declarations outside its measured intervals to
verify that its nonterminal probes remain nonterminal; its whole-process time
must not be compared as a scalar-application speedup.

All frozen source/asset and executable hashes matched before/after the gate.
An independent audit reran the exact catalog migration in four fresh processes,
obtained byte-identical shipped gzip files, and checked adversarial transport
inputs. Migration, load timings and numerical receipts are retained in
`TMP/gamma-output-catalog-rollout.B4yQ1n/`; the independent audit is in
`TMP/compressed-catalog-audit.GINagg/` in the RustRed workspace. Neither these
tests nor the 74-output convention constitute arbitrary-index family closure.

## Historical full-raw-catalog native migration

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

## Weighted vacuum terminal normalization

This normalization milestone used RustRed pin
`8ad62b964de6f3a508fd165dc6ac2509f25839fe`.
The four new sidecars retain the verified positive-power U equivalences and
add exact weighted projections for quadratic numerator terminals. RustRed's
generic service authenticates integer unit-Jacobian momentum symmetries of the
active unit-mass denominators, then uses native Symbolica linear algebra to
express each supported numerator as a scalar/pinch combination modulo those
symmetries. The implementation does not dispatch on these parent names or on
four loops. Strict output descent, declared-terminal binding and one-hop
fixed-point closure are checked before installation.

| Parent | Raw program keys | Prior U outputs | Weighted positive outputs | Projected numerators | Sidecar bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| H | 386 | 52 | 22 | 30 | 8,934 |
| FG | 145 | 26 | 16 | 10 | 5,759 |
| BMW | 179 | 37 | 17 | 20 | 6,291 |
| X | 445 | 64 | 19 | 45 | 9,811 |
| Total | 1,155 | 179 | 74 | 105 | 30,795 |

These are 74 family-local positive outputs, not 74 independent or minimal
masters. At the weighted-normalization milestone all twelve existing
program/catalog/descriptor assets remained byte-identical; the subsequent
value-only pruning above changes catalogs but not programs or sidecars. No
candidate rules or master values were regenerated. Export
and fresh dirty-context reload independently checked every one of the 1,155
raw expansions, all 1,260 output coefficients and both ordered native variable
maps. Each family's retained symmetry generators and affine projection columns
were independently replayed before opening that family's unchanged exact value
catalog (84 generators in total). The corpus has no preparation skips.
Unsupported shapes and exhausted inherited U-preparation bounds can retain
their original keys; exhaustion of the weighted preparation's structural
limits instead returns a typed error.

The sidecar uses a distinct native envelope and is trusted generated data for
this pinned stack. Its decoder regenerates the finite proof and compares every
saved expansion before returning an installable owner; framing alone conveys
no authority. The current catalog check binds the exact final output set;
the normalizer itself still checks every raw program declaration.
The thin Vakint constructor `from_reducer_with_terminal_normalization` installs
that owner into an empty cache. The plain and unit-alias constructors retain
their prior behavior, and no cache is silently cleared. Hot application uses
the core's existing coefficient arithmetic and memoization; Vakint adds no
algebra or reduction engine.

The pinned Vakint gates pass the unchanged 83-case through-three-loop selection,
all fifteen original references, all sixteen expanded-numerator/pinch pairs,
eleven focused constructor/catalog/loader checks and three fixtures. The two
four-loop correctness processes took 20.58/36.04 s wall time, with peak RSS
1,842,588/2,003,816 KiB. They ran sequentially on CPU 83 and include the separate
FMFT oracle; these are not scalar-only timings. The original nonunit mass/scale
inputs, precision and tolerances remain unchanged, with forbidden FORM paths
in both native stages. The matched benchmark below also passes all 54 numerical
comparisons. Focused and inventory filters overlap; these are not additional
distinct physical inputs.

## Historical U-only vacuum terminal normalization

The preceding `2b50267c` runtime pin explicitly steered
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

That production loader opted into `from_reducer_with_terminal_aliases` on a fresh
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
tables are used. The compressed `.rrcat.bin.gz` records are exact PR expressions, **not new
20,000-digit master evaluations**. Requesting more precision than a surviving
master/constant source contains emits the existing warning through Vakint's
logger; it does not increase source accuracy. Missing Laurent orders are errors.

## Offline reproduction

The historical offline commands below use RustRed runtime pin
`8ad62b964de6f3a508fd165dc6ac2509f25839fe`, with the Symbolica license supplied
in the environment. These commands generate fresh programs when intentionally
requested; the shipped migration itself converted the saved programs without
generating new rules. To reproduce the current generated packages, use the
[current producer recipe](../../../rustred_package_refresh.typ) instead.

```sh
cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  generate /path/to/vakint/data/rustred/four_loop/h.toml 9 6 default \
  /path/to/new/h.candidates.rrbin /path/to/new/h.report.toml

cargo run --release --locked --offline --no-default-features \
  -p rustred-app --example candidate_bundle -- \
  verify 10 /path/to/new/h.candidates.rrbin \
  /path/to/offline/raw/h.rrcat.bin

gzip -n -9 -c /path/to/new/h.candidates.rrbin > /path/to/new/h.candidates.rrbin.gz
```

The helper requires new output paths. The zero-based nonpositive index list is
`9` for H/X and `8,9` for FG/BMW; these are auxiliary scalar-product coordinates,
not physical propagators. The worker count above is an explicit example and can
be changed. Generation uses ordinary RustRed IBPs and no FMFT-derived rules.
The fresh-process `verify` operation checks the family/catalog binding and
declared raw terminal keys; it requires a separately prepared full raw diagnostic
catalog, not the shipped output-only catalog. Runtime validation instead binds
the compressed catalog to the installed normalizer. Neither operation is a
symbolic closure certificate.

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

### Matched weighted-sidecar comparison, 20 September 2026

The control is a fresh run of the frozen U-only `2b50267c` executable; the new
runtime pins `8ad62b96` and uses the four independently replayed sidecars. Each
passes the same nine inputs and 54 numerical comparisons, with 108 timed calls.
Both use the preceding public test unchanged, release optimization level 3,
LTO disabled, affinity 88–93, one scalar caller and nested pools capped at one.
The processes run separately, with a 600-second deadline and 8-GiB virtual-address
cap each; actual RSS is measured independently. Other builds/validation use
disjoint CPUs on the same shared host. This is one observed pair of pinned
runtimes, not a statistical speedup distribution or an isolated kernel change.

| Input | First RustRed call (s), U → weighted | Warm RustRed median (ms), U → weighted | Warm FMFT median (ms), U → weighted |
| --- | ---: | ---: | ---: |
| H, first propagator cubed | 6.231 → 5.587 | 84.985 → 80.665 | 1493.497 → 1492.737 |
| H, expanded D7 numerator | 0.021 → 0.021 | 20.887 → 21.106 | 135.081 → 135.994 |
| FG, first propagator cubed | 1.043 → 1.110 | 79.621 → 76.572 | 223.242 → 223.695 |
| FG, expanded D7 numerator | 0.763 → 0.654 | 31.135 → 30.704 | 149.504 → 152.586 |
| BMW, first propagator cubed | 3.594 → 2.941 | 67.392 → 62.268 | 709.840 → 726.586 |
| BMW, expanded D7 numerator | 0.030 → 0.024 | 26.026 → 22.702 | 136.294 → 136.249 |
| X, first propagator cubed | 30.811 → 21.411 | 120.613 → 113.522 | 4780.624 → 5857.101 |
| X, expanded D7 numerator | 0.631 → 0.367 | 52.306 → 53.005 | 163.334 → 196.491 |
| Factorized four-tadpole | 0.015 → 0.016 | 14.751 → 15.309 | 128.703 → 135.580 |

Observed H/X first-parent calls improve by 1.12×/1.44×, including lazy program
load and native sidecar proof reconstruction. This is not a loader-only result.
FG's first cubed-line call regresses, as do the warm H/X expanded-D7 calls and
the factorized control. The substantial movement in the independent FMFT X
warm median also cautions against attributing all timing movement to normalization
on this shared host; no observations were removed.
The dotted public targets have power three, unlike the separate core D1² study.

Warm RustRed medians are below FMFT for all nine tested inputs, but the first
four cubed-parent calls remain slower than FMFT. For H/X the new initial FMFT
calls are 1.479/4.787 s versus RustRed's 5.587/21.411 s. Initial FG/X expanded-D7
calls also remain slower than FMFT. Cached repeats do not predict unseen-target
costs, and 179→74 outputs does not imply a proportional runtime improvement.

| Whole benchmark process | U plan | Weighted sidecars |
| --- | ---: | ---: |
| Wall time | 97.22 s | 89.82 s |
| User + system CPU | 92.97 + 3.63 s | 85.44 + 3.75 s |
| Peak RSS | 2,155,424 KiB | 2,039,496 KiB |

These totals include both backends, initialization and untimed comparisons;
RSS is not a per-backend cache measurement. No rules, raw catalog values or
numerical master precision changed. The frozen executables, full observations,
commands, exact medians and audits are retained in the RustRed workspace's
`TMP/gamma-terminal-normalization-rollout.GN7v1b/`. The earlier comparisons
below remain historical evidence and use different controls.

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

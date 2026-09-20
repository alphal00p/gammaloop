# Four-loop RustRed acceptance and candidate diagnostics

The public `EvaluationMethod::RustRed` now loads the shipped four-loop candidate
programs and terminal catalogs under `data/rustred/four_loop/`. Ordinary
evaluation does not regenerate parametric IBPs or invoke FORM. The programs
retain candidate status: loading and pointwise rule checks do not establish
arbitrary-index closure or constitute a sealed completeness certificate.
Vakint's defaults and existing FORM-backed methods remain unchanged.

Two ignored numerical suites exercise this actual public backend through
`EvaluationOrder::rustred_only()`: the fifteen original numerical references
and sixteen expanded-numerator/pinch pairs. Neither suite injects a test-local
reducer or bypasses public method dispatch. Programs are loaded lazily and
reused; these tests contain no IBP solve. FeynKit tensor reduction and the
RustRed scalar tail run with an invalid FORM path. A separately supplied valid
FORM executable is used only by the independent FMFT comparison peer.

The older `candidate_*` diagnostics remain available for studying explicitly
supplied families and generation timings. They are distinct from the public
acceptance suites and deliberately regenerate their candidate rules.

## Fifteen original numerical references

`public_acceptance::all_fifteen_numerical_references_match_fmft_and_public_rustred`
compares FMFT and public RustRed with all fourteen original analytic four-loop
references and the decorated four-tadpole case whose historical name contains
`1l`. It preserves the original inputs, mass and external-momentum values,
spectator symbols, normalization, five Laurent terms, 32-digit arithmetic,
relative reference tolerance `1e-30`, and maximum pull `1.0`. The fixture-only
tests separately check inventory and literal-input agreement with the originals.

```sh
export SYMBOLICA_LICENSE="<your current Symbolica license>"
export VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH=/absolute/path/to/form
cargo test -p vakint --release --locked \
  --test rustred_four_loop_tests \
  all_fifteen_numerical_references_match_fmft_and_public_rustred \
  -- --ignored --nocapture --test-threads=1
```

Optionally set `VAKINT_4L_CANDIDATE_FAMILY_FILTER` to `H`, `FG`, `BMW`, or `X`.
Selection follows the retained matcher parent, not the historical test name.
The complete public inventory passed on 2026-09-19: **15/15 references**,
132.59 s full-process wall time (132.06 s test body). The public integration
target does not require `experimental-rustred`. These timings include program
loading, both reduction peers and numerical checks; they are not scalar-only
timings or a RustRed-versus-FMFT speed comparison.

The separate seven-test default-backend/routing regression selection also
passed (1.38 s test body), covering opt-in dispatch, mass/power admission,
retained K6 routing, negative orientation and native master normalization.

## Runtime-generation diagnostics

The ignored `candidate_parent_dotted_and_pinch_match_fmft` test has three phases:

1. RustRed solves the supplied family and applies candidate rules to a scalar
   parent, a dotted parent, and a pinch. RustRed owns exact point guards,
   descent, memoization, and unresolved-target errors. At least one actual rule
   application is required; a terminal-only comparison is labeled separately.
2. An explicitly enabled offline FMFT oracle evaluates the declared fixed
   residuals. These values form a frozen terminal catalog.
   They are not used as IBP rules or as additional master declarations.
3. The existing comparative harness runs FeynKit, the legacy FMFT peer, and the
   candidate RustRed scalar peer. The candidate peer has an invalid FORM path
   and uses Vakint's existing pure-Rust FMFT master finalizer.

Tensor/scalar-product terms are routed through the RustRed-owned
`FamilyScalarNumeratorService` before candidate rule application.  The adapter
does not expand scalar products itself: each lowered key, exact spectator, and
common-mass power is carried into the same candidate reducer.  This exercises
the FORM-free scalar tail after the existing FeynKit tensor prepass, while the
finite-target scope remains unchanged.

The default CSV under `tests/inputs/experimental_four_loop_h.csv` is test input,
not an engine dispatch. Its rows define physical graph edges and momentum
coordinates; the auxiliary slot has edge `(0,0)`. Both the RustRed family and
oracle input are constructed from the same supplied descriptor. Another input
may be supplied with `VAKINT_4L_CANDIDATE_PARENT_INPUT`.
The adjacent FG, BMW and X CSV inputs use the same descriptor schema, so the
same compiled test binary can study them without recompilation.
`VAKINT_4L_CANDIDATE_DOT_POWER` changes the dotted target's first physical
power (an integer at least two; default two). This is target input, not solver
dispatch. For example, if all three default probes are already declared
terminals, a larger explicit power can exercise an actual recurrence without
weakening the mandatory rule-application assertion.

To run this separate generation diagnostic:

```sh
export VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH=/absolute/path/to/form
export VAKINT_4L_CANDIDATE_WORKERS=1
cargo test -p vakint --release --locked --features experimental-rustred \
  --test experimental_rustred_4l_tests \
  candidate_parent_dotted_and_pinch_match_fmft \
  -- --ignored --nocapture --test-threads=1
```

Supply Symbolica's license through the environment as usual. No license is
stored in this fixture. The explicit FORM executable is used only by the
offline preparation and legacy comparison peer, not by candidate application.

The test prints search time, fixed-residual count, per-target rule counts,
offline terminal expressions, and finite-target parity results. An optional
`VAKINT_ACCEPTANCE_EXACT_DIAGNOSTICS=1` pass also inspects exact Laurent
coefficients before numerical master substitution, without changing the
compared values or tolerances. Catalogs containing floating coefficients are
rejected before exact PR normalization.

After that unchanged numerical gate succeeds,
`VAKINT_4L_CANDIDATE_BENCH_REPEATS=3` additionally measures identical scalar-tail
and full FeynKit-plus-scalar evaluations. One-off candidate search and offline
terminal preparation are excluded. Native memoization is cleared outside each
logical-cold interval; immediate warm-cache observations are reported separately
with their rule counts, not presented as fresh reduction. Lane order alternates
between repeats, and every timed output is numerically checked outside the
timer. FORM process startup and temporary I/O remain part of FMFT evaluation.
Linux CPU deltas include waited-for child processes; RSS fields are explicitly
parent-process before/after snapshots, not phase peaks or FORM-child memory.

Historical probe results from 2026-09-16 are retained here for provenance:
the H parent, dotted parent and pinch passed the unchanged strict comparative
harness, with 26,956 invalid-FORM scalar-tail applications after cache clear;
the same binary passed the three FG targets with 3,362 applications and the
three X targets with 82,637 applications.
These experiments use exact offline projections of only the actually reached
fixed residuals, not guessed master declarations or oracle-supplied IBP rules.
The historical BMW dotted power-three probe passed with 16,113 applications.
Together these are twelve historical parent/dotted/pinch comparisons, not a
measurement of the newer fifteen-entrypoint analytic inventory.

The ignored `candidate_all_four_loop_parents_match_fmft` test reuses this
same comparative harness for all four checked-in descriptors (`H`, `FG`,
`BMW`, and `X`), with parent, dotted, and pinch probes for each descriptor.
Each descriptor gets a fresh RustRed candidate reducer and an offline catalog
containing every finite terminal declared by that run. A matrix pass demonstrates
routing, FeynKit preprocessing, candidate application, and numerical parity
for those finite probes only; it does not certify arbitrary-index closure. It
does not exercise loading through the public method; use the public suites
above and below for that acceptance boundary.

The harness ships complete finite catalogs for the four checked-in parent
descriptors under `data/rustred/four_loop/`; these are
loaded by default and are bound to each RustRed family fingerprint and arity.
Set `VAKINT_4L_CANDIDATE_CATALOG_DIR` to override that bundle and persist a
fresh exact offline value set. The harness writes `h.rrcat.bin`, `fg.rrcat.bin`,
`bmw.rrcat.bin`, or `x.rrcat.bin` atomically after the FMFT phase; each file contains
every finite terminal declared by that candidate and is marked complete. With
`VAKINT_4L_CANDIDATE_CATALOG_ONLY=1`, an existing catalog is required and the
test uses an invalid FORM path for the RustRed/FeynKit scalar tail. This is an
experimental value-cache check, not a completeness proof or IBP artifact.
The same canonical catalog directory serves the public loader; there is no
separately maintained test-catalog copy.
Catalogs produced by the earlier probe-only harness are marked partial and
are intentionally rejected by this strict complete-terminal reload path.
The native catalog codec lives in RustRed and uses Symbolica's Atom/state
serialization; the former line-oriented schema is no longer a runtime input.
The four shipped files preserve all 1,155 exact values and keys and total
30,307 bytes. Loading accepts trusted generated native data, not hostile files.
The `d51721b6` native-catalog runtime gate passes the fifteen-reference and
sixteen-pinch public suites with the existing tolerances and forbidden FORM
paths in the native stages. Retained evidence is under
`TMP/gamma-native-catalog-migration.jP7cXc/` in the RustRed workspace, including
the independent `public-audit.md`. No alias optimization was activated for this
format-only comparison.

The later `f91c47ab` runtime explicitly installs RustRed's verified terminal
routing plan only in the production four-loop loader. The old experimental
`NativeCandidate::from_reducer` keeps raw outputs and accepts populated caches;
`from_reducer_with_terminal_aliases` is the explicit fresh-owner constructor.
All raw catalog declarations and saved program bytes are unchanged. The full
fifteen-reference and sixteen-pinch gates pass again with forbidden FORM paths
in the native stages, at 24.71/47.88 s whole-process wall time. Nine focused
checks include both constructor conventions, and the unchanged 83-test lower-
loop selection also passes. The same nine-input public timing test passes all
54 numerical comparisons; its full table and first-use/warm-cache boundaries
are in the [asset README](../../data/rustred/four_loop/README.md#public-scalar-timing-probe).
Evidence, exact binary hashes, commands and independent audit are retained at
`TMP/gamma-terminal-alias-rollout.4Zi7ej/` in the RustRed workspace.

The subsequent `2b50267c` pin changes only that explicit fresh-owner factory to
RustRed's native U-polynomial equivalence plan with default preparation bounds.
The raw-key binding, plain constructor, public options, assets and numerical
assertions are unchanged. All 83 lower-loop cases, fifteen original four-loop
references, sixteen expanded-numerator/pinch pairs, nine focused checks and
three fixtures pass again. The independent public four-loop gates took
22.56/46.90 s whole-process wall time and retained forbidden FORM paths in both
native stages. These correctness observations include the separate FMFT oracle.
The frozen routing-only executable and the U-enabled executable each pass the
same nine-input timing test, with all 54 numerical comparisons. The measured
table, warm-cache and first-use boundaries are in the asset README; retained
commands, logs, hashes and audit are under
`TMP/gamma-terminal-u-rollout.oFoQwr/` in the RustRed workspace. No rules or
catalog values were regenerated.

The current `8ad62b96` pin loads four additional `.rrnorm.bin` sidecars through
the explicit `from_reducer_with_terminal_normalization` constructor. The core
rebuilds the exact finite proof on cold loading and installs weighted terminal
rows before memoization; Vakint implements no new reduction arithmetic. The
plain and U-alias constructors, original raw catalogs, saved candidate programs,
default evaluation order and numerical assertions remain unchanged. Exact
export and fresh dirty-context reload checked all 1,155 raw expansions and all
105 supported numerator projections, giving 74 family-local positive outputs
(H 22, FG 16, BMW 17, X 19), not a minimal-master or closure certificate.
The pinned gate passes all 83 lower-loop selections, all fifteen original
references, all sixteen expanded-propagator/pinch pairs, eleven focused checks
and three fixtures. The two public correctness processes take 20.58/36.04 s
wall time including FMFT, with unchanged nonunit scales, tolerances and invalid
FORM paths in the native stages. A fresh U-only control and weighted runtime
each pass all 54 comparisons in the same nine-input benchmark. Whole-process
wall times are 97.22/89.82 s and include both backends; the asset README retains
every first/warm observation, regressions and shared-host caveats. Evidence is retained at
`TMP/gamma-terminal-normalization-rollout.GN7v1b/` in the RustRed workspace;
the [asset README](../../data/rustred/four_loop/README.md#weighted-vacuum-terminal-normalization)
describes the native proof and unchanged raw-catalog boundary.

The existing FMFT finalizer expands exact Laurent coefficients before
approximate table substitution. This prevents exact cancelling coefficients
from becoming tiny floating poles; no zero tolerance, epsilon depth or source
table is changed. The historical nine-check finalizer unit run passed,
including genuine missing-order rejection and exact-cancellation regression.

For bounded matrix runs, set `VAKINT_4L_CANDIDATE_FAMILY_FILTER` to one of
`H`, `FG`, `BMW`, or `X` to run only that descriptor. The default remains the
complete four-family matrix; the filter is intended to generate or replay one
catalog at a time without changing the candidate or catalog semantics.

## Historical diagnostic measurements

These measurements predate the current public-backend matrix and the FG
signed-momentum correction described below. They are retained for provenance,
not substituted for fresh acceptance results.

The `candidate_h_rank_four_numerator_matches_fmft` lane uses the same
H tensor numerator and external momenta as the analytic four-loop test, with
`muvsq=3` and `mursq=5`. Its regenerated catalogue contains 386 declared
terminals. A fresh strict catalog-only replay loaded all 386 terms, used an
invalid FORM path, passed parent/dotted/pinch/rank-four probes, and applied
26,956 rules in 441.43 s. The ordinary FMFT comparison in that historical
rank-four lane was stopped after 10m22 without a parity result. No newer
ordinary rank-four H parity result is claimed here. The
`candidate_fg_clover_numerator_case_matches_fmft`
lane reuses the
registered FG parent witness for all four literal Clover probes (unit mass,
non-unit masses, dotted `k1^2`, and the rank-four numerator). The scalar
probes passed finite-target parity in that historical run. The rank-four probe
uses the exact
upstream `user_space::A/B/C` numerator and the dotted first propagator; the
adapter restores routed indexed components to scalar products before the
RustRed-owned numerator lowerer and rejects any residual loop-vector
component. The release run on 2026-09-18 completed in 66.06 s with 3,366 cold
rule applications and an invalid scalar FORM path. A standalone four-slot
clover descriptor is intentionally not shipped: the current matcher witness
represents that case as a contraction of an eight-slot parent.

Historical release-binary parent runs on 2026-09-18 (six candidate workers)
produced the bundled catalogs
with the following measurements. They are retained for comparison only:

| parent | fixed residuals | candidate search | total candidate/oracle/parity run | dotted applications |
| --- | ---: | ---: | ---: | ---: |
| FG | 145 | 158.23 s | 278.17 s | 3,362 |
| BMW | 179 | 317.98 s | 472.24 s | 16,113 |
| X | 445 | 145.31 s | 430.49 s | 82,637 |

The historical H rank-four catalog-only replay took 441.43 s with 386 terms and
26,956 applications. A bundled FG catalog-only replay with an invalid FORM
path completed in 50.38 s (3,362 applications). These timings include the
finite candidate/oracle harness boundary and are not claims about a complete
four-loop artifact or production RustRed evaluation.

The older `candidate_all_four_loop_analytic_acceptance_cases_match_fmft`
generation diagnostic contains the same **fifteen numerical acceptance
entrypoints**, spanning the
H, FG, BMW and X descriptors, registered numerator and non-unit-mass cases,
and the historically mislabelled four-tadpole entrypoint. Inventory
implementation and input validation do not imply that these numerical runs
have all passed. A previous serial filtered-results table has been withdrawn:
its complete run logs were not available for independent verification. The
historical twelve-probe evidence and the explicitly limited H rank-four
catalog-only evidence above remain separate from this new inventory.

The inherited environment of a later six-worker attempt lacked
`SYMBOLICA_LICENSE`; Symbolica consequently entered restricted mode and
reported `requested 6 sector workers, but Symbolica permits 1 on execution
threads`. This is not evidence that the user's supplied license is limited to
one worker. Export the current license before running either the serial or
parallel matrix; no new six-worker timing is inferred from that failed setup.

## Expanded propagator numerators versus explicit pinches

The public `propagator_pinches.rs` lane defines sixteen paired-input checks: four
for each of H, FG, BMW and X. Numbering propagators from one, the expanded
numerator is respectively `D1`, `D7`, `D1*D2`, or `D5*D6`, where each `Di`
includes the mass term with the input propagator's sign convention. The
paired input has no such numerator and instead contracts the corresponding
one or two propagators. Scalar products and products of denominators are
expanded, so the test exercises numerator handling, not a textual cancellation
of identical factors. The selected double pinches preserve four-loop,
non-scaleless cases.

Both inputs are matched independently through Vakint's existing topology
matcher; their retained routing witnesses select the shipped program. A
pinch may therefore use a different registered parent from the numerator
input. The harness uses the FeynKit tensor prepass and tests both the FMFT peer
and the public RustRed scalar peer, including comparison between the
peers. Both FeynKit and the RustRed peer use an invalid FORM path. These checks
use
`muvsq=3`, `mursq=7`, 32-digit numerical arithmetic, and relative tolerance
`1e-20`, making missing mass terms or powers observable.

Run the numerical matrix explicitly:

```sh
export SYMBOLICA_LICENSE="<your current Symbolica license>"
export VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH=/absolute/path/to/form
cargo test -p vakint --release --locked \
  --test rustred_four_loop_tests \
  sixteen_numerator_propagator_pinches_match_fmft_and_rustred \
  -- --ignored --nocapture --test-threads=1
```

Set `VAKINT_4L_CANDIDATE_FAMILY_FILTER=H`, `FG`, `BMW`, or `X` to select the
four numerator inputs belonging to one parent. An explicitly pinched partner
can still load a program for another matched parent. Public programs and
terminal catalogs are reused; neither IBP generation nor offline catalog
generation occurs in this suite. On 2026-09-19 the complete public matrix
passed **16/16 pairs**, checking each identity in both FMFT and RustRed and
checking cross-backend agreement. Full-process wall time was 208.69 s
(208.10 s test body). This includes both peers and program loading, not just
RustRed scalar application. The earlier runtime-generation version's failure
is retained below as historical evidence. The unignored
`sixteen_inputs_are_expanded_and_have_registered_four_loop_witnesses`
test checks fixture construction separately and cannot establish numerical
parity by itself.

Historical checkpoint on 2026-09-18: the expanded-input/matcher fixture test
passed for all sixteen pairs (2.42 s). The six-worker numerical run completed with a
failure after 67.37 s: all four H pairs and FG/D1 passed both the within-backend
and cross-backend checks. FG/D7 versus its single pinch then failed **inside
FMFT**, before the RustRed lane: the real epsilon^-4 coefficients were 1.5 and
0. The other ten pairs were not reached. The H double pinch D5*D6 exercised
12,554 new RustRed recurrence applications; the other passing pairs used
declared terminals. No tolerance was weakened. Investigation after resuming
identified a signed-momentum routing error in the FMFT FG adapter: the fourth
loop momentum must map to `-p6`, not `p6`. The adapter correction and an exact
routing regression are implemented. The corrected public matrix passed as
reported above, after refreshing the four affected FG catalog entries. This
historical failure was not a RustRed reduction failure.

The passing sixteen-pair matrix establishes these concrete cancellation and
routing identities, not arbitrary-index family closure or comprehensive
coverage of every four-loop input.

## Offline shipped-catalog audit

`catalog_validation::shipped_terminal_projections_match_fmft` evaluates only
the already-declared terminal keys with FMFT and compares their exact master
projections with `data/rustred/four_loop/*.rrcat.bin`. It does not solve IBPs or
declare additional terminals. With the license and explicit oracle environment
set as above, run this separate offline diagnostic:

```sh
cargo test -p vakint --release --locked --features experimental-rustred \
  --test experimental_rustred_4l_tests \
  shipped_terminal_projections_match_fmft \
  -- --ignored --nocapture --test-threads=1
```

The family filter is supported. On 2026-09-19, complete re-evaluation of all
**1,155 terminal keys** finished in 199.29 s. H (386), BMW (179), and X (445)
had zero differences; FG (145) had four projections changed by the signed
routing correction. Those four entries were refreshed, and a subsequent FG
audit passed with **zero differences** in 23.62 s. No parametric IBP program
was regenerated for this catalog correction.

Without `VAKINT_4L_CATALOG_AUDIT_OUTPUT`, every projection difference fails the
audit. Supplying that variable requests an explicit offline refresh into the
chosen directory instead; such a refresh is not a passing validation of the
previous catalogs. Review updated values, rebuild embedded assets, and rerun
the public numerical suites before reporting parity.

Historical diagnostics and the new public tests have finite target scope.
Neither proves arbitrary-index closure, supplies a sealed four-loop artifact,
nor establishes twenty-thousand-digit master accuracy. Shipping candidate
programs removes runtime rule generation; it does not change those proof and
precision limits.

# Experimental four-loop candidate comparison

This opt-in development lane can compare finite integral targets before the
candidate rule family has passed RustRed's complete artifact certification.
It does **not** register a four-loop production artifact, modify Vakint defaults,
or establish arbitrary-index closure. Production `EvaluationMethod::RustRed`
continues to require shipped certified artifacts.

The ignored test `candidate_parent_dotted_and_pinch_match_fmft` has three phases:

1. RustRed solves the supplied family and applies candidate rules to a scalar
   parent, a dotted parent, and a pinch. RustRed owns exact point guards,
   descent, memoization, and unresolved-target errors. At least one actual rule
   application is required; a terminal-only comparison is labeled separately.
2. An explicitly enabled offline FMFT oracle evaluates only the declared fixed
   residuals actually reached. These values form a frozen terminal catalog.
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

Once the local RustRed candidate bridge is available, run:

```sh
export VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH=/absolute/path/to/form
export VAKINT_4L_CANDIDATE_WORKERS=1
cargo test -p vakint --release --features experimental-rustred \
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

On 2026-09-16 the H parent, dotted parent and pinch passed the unchanged strict
comparative harness. The invalid-FORM scalar tail performed 26,956 new rule
applications after its cache was cleared; parent and pinch themselves are
declared terminals. A separate source/evidence audit confirmed the pass. The
same binary also passed the three FG targets, with 3,362 new applications.
The three X targets also passed, with 82,637 new applications.
These experiments use exact offline projections of only the actually reached
fixed residuals, not guessed master declarations or oracle-supplied IBP rules.
The BMW parent, dotted power-three target and pinch now use the same harness;
the dotted target passed with 16,113 new applications.  The power-two BMW
probe is a declared residual and is therefore not counted as recurrence
evidence.

The ignored `candidate_all_four_loop_parents_match_fmft` test reuses this
same comparative harness for all four checked-in descriptors (`H`, `FG`,
`BMW`, and `X`), with parent, dotted, and pinch probes for each descriptor.
Each descriptor gets a fresh RustRed candidate reducer and an offline catalog
containing every finite terminal declared by that run. A matrix pass demonstrates
routing, FeynKit preprocessing, candidate application, and numerical parity
for those finite probes only; it does not certify arbitrary-index closure or
enable four-loop `EvaluationMethod::RustRed`. The production registry remains
sealed-artifact-only until a complete authenticated four-loop artifact and
terminal manifest are available.

Set `VAKINT_4L_CANDIDATE_CATALOG_DIR` to persist exact offline values for later
FORM-free candidate reruns. The harness writes `h.rrcat`, `fg.rrcat`,
`bmw.rrcat`, or `x.rrcat` atomically after the FMFT phase; each file contains
every finite terminal declared by that candidate, is marked complete, and is
bound to the RustRed family fingerprint and arity. With
`VAKINT_4L_CANDIDATE_CATALOG_ONLY=1`, an existing catalog is required and the
test uses an invalid FORM path for the RustRed/FeynKit scalar tail. This is an
experimental value-cache check, not a completeness proof, IBP artifact, or
production master catalog.
Catalogs produced by the earlier probe-only harness are marked partial and
are intentionally rejected by this strict complete-terminal reload path.

The existing FMFT finalizer now expands exact Laurent coefficients before
approximate table substitution. This prevents exact cancelling coefficients
from becoming tiny floating poles; no zero tolerance, epsilon depth or source
table is changed. All nine finalizer unit checks pass, including genuine
missing-order rejection and the new exact-cancellation regression.

For bounded matrix runs, set `VAKINT_4L_CANDIDATE_FAMILY_FILTER` to one of
`H`, `FG`, `BMW`, or `X` to run only that descriptor. The default remains the
complete four-family matrix; the filter is intended to generate or replay one
catalog at a time without changing the candidate or catalog semantics.

The added `candidate_h_rank_four_numerator_matches_fmft` lane uses the same
H tensor numerator and external momenta as the analytic four-loop test, with
`muvsq=3` and `mursq=5`. Its regenerated catalogue contains 386 declared
terminals. A fresh strict catalog-only replay loaded all 386 terms, used an
invalid FORM path, passed parent/dotted/pinch/rank-four probes, and applied
26,956 rules in 441.43 s. The ordinary FMFT comparison was stopped after
10m22 without a parity result, so no ordinary rank-four parity claim is made
yet. The new `candidate_fg_clover_numerator_case_matches_fmft` lane reuses the
registered FG parent witness for all four literal Clover probes (unit mass,
non-unit masses, dotted `k1^2`, and the rank-four numerator). The scalar
probes pass finite-target parity. The rank-four probe now uses the exact
upstream `user_space::A/B/C` numerator and the dotted first propagator; the
current RustRed/FeynKit tail reaches the expected lowered expression but
fails closed because the evaluator has no parameter-map entry for the
remaining tensor scalar `dot(k(3),k(3))`. This is a pending tensor-bridge
issue, not a catalog or closure failure: the optimized run on 2026-09-18
completed the six scalar probes in 82.16 s and exited 101 on that final
probe, so no rank-four numerical parity claim is made. A standalone four-slot
clover descriptor is intentionally not shipped: the current matcher witness
represents that case as a contraction of an eight-slot parent.

This checkpoint is **twelve finite-target comparisons plus one strict H
rank-four catalog-only replay**, not all fifteen existing
four-loop numerical acceptance entrypoints or the nineteen registered graph
classes. It does not prove arbitrary-index closure, ship a four-loop production
catalog, or establish twenty-thousand-digit master accuracy. Tensor-bearing
candidate acceptance and non-unit-mass candidate comparisons remain separate
gates.

= Compact energies and factor collection in soft counterterms

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

The square-root function is now named `OnShellEnergy(owner, s) = sqrt(s)`;
the recorded runs called it `Esurface`. The prose uses the current name except
where explicitly quoting the original component-argument syntax. Today,
`Esurface(a1, ..., an) = a1 + ... + an` denotes the surrounding signed energy
sum, including energy shifts. This naming change does not alter the historical
measurements below.

#set text(size: 10pt)

#strong[September 25 follow-up:] Evaluator construction now restores main's
Spenso scalar aliases: scalar bodies of at least 4096 bytes are eligible, with
branch and orientation guards kept visible. Aliases are introduced after network
parsing, at the evaluator boundary; Taylor-stage compact energies retain their
explicit dependencies. The stack also includes the batched scalar-sum parser fix
and the newer Symbolica series implementation. The measurements and descriptions
below record the original September 24 configuration, with automatic aliases
disabled; they are not timings of this follow-up.

#strong[Single-invariant follow-up:] The current energy call is
`OnShellEnergy(owner, s)`, with `s = m^2 + qx^2 + qy^2 + qz^2` and
`dE/ds = 1/(2E)`. Routed momenta, mass scopes, and all Taylor parameters remain
explicit inside `s`; only one function argument participates in the chain rule.
The owner is retained until final evaluator normalization. The hard/soft
rescaling now divides the whole invariant by `r^2`, with the same UV mass
completion as before. The evaluator defines the function as `sqrt(s)`.
Canonicalization flips signs inside squared scalar factors and distributes
numeric coefficients only. It does not expand squared momentum sums or graph
numerators. Taylor probe orders, factor collection, and evaluator settings are
unchanged in this comparison. The representation and timings below describe the
original component-argument implementation.

For the bubble-dressed double triangle in localized direct 3D, a matched bounded
run stopped after the same seven Taylor operations. Total wall time fell from
397.9995 s to 26.6402 s, and peak process RSS from 9.1793 GiB to 0.1057 GiB.
The worst numerator's absolute-order-zero probe still produced all seven
coefficients, but fell from 314.6827 s to 0.01369 s; its coefficient payload fell
from 6,513,098 to 167,561 bytes. The three coefficients actually retained took
0.00282 s. The 102 root probes and root expansions in that operation now account
for 6.2683 s and 8.0222 s respectively. These are prefix measurements, not a
completed diagram-generation timing. Frozen binaries, source hashes, commands,
resource samples, structured logs and comparison data are under
`/tmp/soft-ct-invariant-energy-2026-09-25/`; the baseline is
`/tmp/soft-ct-expression-probe-2026-09-25/`.

The follow-up passed 36 focused tests and two additional direct/projected route
oracles. The existing depth-three nested-banana discrepancy was reproduced with
the same values as the prior recorded run. Quark self-energy generation and UV
profiles passed in localized 3D, explicit 3D, and projected 4D. Clippy passed with
the existing `uv/profile.rs` argument-count warning. A separate full-generation
attempt, concurrent with validation builds, stopped at its 180-second limit:
ten Taylor operations completed, with peak RSS 301.1 MiB. Operation eleven's
92 root probes took 55.92 s; its root expansions were still in progress. This
attempt did not reach evaluator construction. The adjacent
`soft-ct-compact-energies-2026-09-24/invariant-comparison-2026-09-25.json`
retains the run identities, measurements and validation summary.

#strong[Probe-reuse follow-up:] The next change retains the branch-probe
coefficients, then composes normalized numerator polynomials into them. Energy
denominators are differentiated only during the probe. Ten of eleven numerator
probes also supplied the required coefficients without another expansion.
The matched first seven operations took 19.6314 s instead of 26.6402 s.
Operation seven's final branch-series work fell from 8.0222 s to 1.9346 s;
its probe time remained comparable at 6.4907 s versus 6.2683 s.
Peak RSS rose from 108.3 MiB to 115.8 MiB because probes remain stored until used.

Thirty-seven focused tests and two additional route checks passed. The new
regression verifies the denominator derivative count as well as the exact
Taylor result. Nested-top generation and UV profiles passed in all three
routes: 21 summed limits per route and 294 orientation checks in localized 3D.
The known nested-banana discrepancy is unchanged. The separate 180-second
full-generation attempt completed ten operations and 26 of operation eleven's
92 final branch compositions, with peak RSS 546.9 MiB. All 92 probes in that
operation completed and took 62.41 s. This attempt ran concurrently with Clippy
and bounded nested-top validation; it is not a controlled end-to-end timing.
The diagram remains unfinished. Measurements and identities are retained in
`soft-ct-compact-energies-2026-09-24/reused-root-probes-2026-09-25.json`;
full artifacts are in `/tmp/soft-ct-reused-root-probes-2026-09-25/`.

#strong[September 26 Symbolica Horner fix:] Symbolica, Numerica and Graphica
are pinned together to main revision
`70375b9ef06411d8bac3e14d1f57382384547eef`. Its Horner implementation groups terms
without retaining successive copies of the remaining sum. The update also
requires `set_precision` implementations: complex and scalar wrappers delegate
to their components; the const-precision `VarFloat` keeps its precision fixed.
The local factor collector, opaque function/power policy, scalar alias policy,
and evaluator settings are unchanged.

The same scalar-only Horner reproducer preserves exact results. For 250, 500 and
1,000 factor groups, Horner takes 1.107, 2.501 and 6.256 milliseconds, compared
with 43.623, 164.283 and 619.059 milliseconds before the update. At 1,000 groups,
the Horner checkpoint's peak RSS falls from 797.5 MiB to 19.0 MiB.

The saved 229-MB evaluator input now completes its isolated collector replay.
Horner takes 3.849 s, factor collection 5.131 s, and the full process including
parsing takes 40.647 s with peak sampled RSS 1.960 GiB. The resulting atom is
227,579,143 bytes; the original reparsed atom is 229,126,703 bytes. The old Horner
pass exceeded the 5-GiB diagnostic limit before returning. These are scalar
collector measurements, not complete diagram timings. Raw logs and reproducer
sources are retained in `/tmp/soft-ct-symbolica-horner-fix-2026-09-26/`.

The localized-3D Double triangle + bubble retry uses the same 900-second and
12-GiB limits, with no concurrent task build. All 33 Taylor operations finish
by 494.44 s; integrand assembly finishes at 552.13 s. Evaluator preprocessing
then completes at 889.90 s, including the previously failing final scalar
collection (30.52 s from network execution to alias retention). The resulting
scalar is 227,635,485 bytes. The process reaches evaluator expression
preparation and is stopped by the wall limit at 900.39 s, with peak sampled RSS
3.575 GiB. The preceding revision exceeded 12 GiB before completing scalar
collection, reaching a sampled peak of 12.592 GiB. Generation still does not
produce a completed evaluator or UV profile within the cap. These are single
runs on a shared host, not a controlled whole-pipeline speedup comparison.

The CLI builds and 92 focused tests pass (85 compact-energy, series, direct-3D,
final-integrand and evaluator checks, plus seven numeric-wrapper tests).
Clippy passes with its existing `too_many_arguments` warning in `uv/profile.rs`;
formatting, workspace-hack verification and CI metadata generation pass.
Source/binary identities and measurements are recorded in
`soft-ct-compact-energies-2026-09-24/symbolica-horner-fix-2026-09-26.json`.

#strong[September 26 guarded scalar sharing:] Numerator-family lowering now
uses direct function-head/tag lookup and simultaneous exact formal binding.
Duplicate formal arguments are rejected before preparation. After outer tensor
contraction and factor collection, repeated complete scalar functions and powers
are shared through `AliasedAtom` when definitions plus references save storage.
Extraction stops at each whole function/power and keeps its body unchanged.
Exact parameters, existing aliases, and branch guards remain visible. Candidates
containing algebraic parameter or function-tag keys are excluded because
evaluator Hornering can change their exact identity inside aliases. Recorded
ordinary function arguments, including `OnShellEnergy` invariants, remain eligible.
This operation occurs after Taylor expansion and formal-argument pruning.

The first prototype, using recursive native alias extraction, completed localized
3D double-triangle-plus-bubble generation in 824.722 s at 3.314 GiB peak RSS;
the previous same-Symbolica run reached its 900 s cap before evaluator building.
Direct family-call lowering fell from 121.850 s to 2.358 s. That prototype reduced
the scalar from 227,635,485 to 39,948,191 bytes and built the Symbolica evaluator
in 47.271 s. These measurements belong to the initial prototype, before the
exact-key safeguards; they must not be attributed to the final extraction rule.

A scalar-only probe, before applying the evaluator's exact parameter/tag catalog,
reduces 227,579,143 bytes after Horner/factor collection to 55,107,780 bytes
including 143 definitions, with exact reconstruction in one replacement pass. Extraction
takes 1.397 s. A blanket exclusion of all algebraic function arguments loses
most of this saving: every excluded candidate contains an `OnShellEnergy` invariant,
which is an ordinary function argument. Parameter keys and fixed function tags
therefore require distinct handling. The first prototype's full UV profile
reached its separate 900 s limit without producing a completed result; numerical
UV certification is not established by successful evaluator generation.

The final implementation passes all 89 focused tests, formatting, Clippy, and the
CLI build. Clippy retains the pre-existing `uv/profile.rs` argument-count warning.
Its fresh localized-3D run uses the same settings and 900 s / 12 GiB limits, with
no concurrent task build. It creates 1,247 aliases, reducing 227,635,485 bytes to
49,684,353 bytes in 1.638 s. Direct family lowering takes 2.387 s, and total
evaluator preprocessing takes 228.539 s versus the baseline's 337.736 s.
Family combination/preparation now accounts for 162.373 s of this total.

This final run reaches its 900.895 s wall cap during Symbolica construction,
with peak RSS 3.426 GiB. Integrand assembly took 629.408 s, versus 541.846 s in
the first prototype run. Selector/expression preparation then took 20.249 s;
Symbolica building started at 878.291 s with a 55,236,632-byte aliased input.
The final run therefore does not establish a completed end-to-end generation
time or full-graph numerical result. The earlier prototype's orientation-zero
and 102-orientation summed smoke checks both returned finite nonzero values,
but these do not replace final-source numerical or UV validation.
Source/binary identities, stage timings, caps, and validation results are in
`soft-ct-compact-energies-2026-09-24/guarded-scalar-sharing-2026-09-26.json`.

Comparing the prototype and final 3D runs, all recorded Taylor payload streams
match exactly after excluding timing and source-location metadata: expressions,
byte sizes, orders, scopes, row order, and counts. Both execute 33 Taylor
operations, 2,201 root probes, 1,984 root series, 223 numerator-body probes/series,
and 2,424 prepared rows. The extra assembly time is distributed across Taylor
calls (+37.021 s), input preparation (+11.070 s), inter-operation assembly gaps
(+28.942 s), and final assembly (+10.339 s). The logs demonstrate no observed
symbolic-work growth before evaluator preparation; they do not record sufficient
CPU or host-load information to attribute the timing variation to contention.

The same frozen final binary was then run on the projected-4D route: 4D Taylor
counterterms, final 3D integrand, explicit orientation sum. Assembly took
642.388 s; evaluator preprocessing took 74.701 s, including 0.825 s for 47,202
family calls. Expression preparation took 1.305 s, so Symbolica construction
started at 718.970 s with a 66,255,467-byte input. The process crossed its 12 GiB
RSS limit and was stopped at 790.598 s, with sampled peak 12.513 GiB. This is a
memory-limited builder run, unlike the time-limited localized-3D run.

This run exposes a coverage gap in scalar sharing: Spenso first stores the
entire 66,255,441-byte scalar in one alias, leaving an 8-byte root. The new pass
walks the root and retains old alias bodies unchanged, so it creates zero new
aliases and does not compact this large definition. Earlier root factor
collection also sees only the handle. Extending sharing into
existing definitions while preserving exact keys requires a further change;
whether it resolves builder memory growth remains untested.

A bounded 2.84 s userspace stack sample during the earlier 4D construction phase
recorded 1,243 samples with no losses. 79.28% of inclusive sampled cycles pass
through `build_projected` and generic `replace_multiple`; 45.76% and 39.59% of
self cycles are in `Pattern::could_match` and `AtomView::replace_no_norm`.
These overlapping percentages identify an earlier substitution hotspot, not the
cause of the later evaluator memory growth. The sample does not distinguish the
individual substitution callsites inside `build_projected`. No process remained
running and no completed 4D evaluator or numerical/UV profile was produced.

#strong[September 27 retained-definition compaction:] A new jj change,
`mvpltvut`, closes the coverage gap above. The existing opaque factor collector
now visits each retained Spenso scalar body before sharing. Sharing counts the
root and every stored definition once, rewrites their complete function/power
factors together, and preserves existing keys. New factor definitions remain
verbatim; the alias dependency graph is not unfolded. Evaluator settings and
Symbolica revision `70375b9ef06411d8bac3e14d1f57382384547eef` are unchanged.

All 91 focused tests pass, including a large scalar hidden in a Spenso definition
and a nested alias graph with shared singular factors. The latter checks exact
reconstruction, a stable second pass, and inactive singular branches in eager and
SymJIT evaluators. Formatting, Clippy and the CLI build pass; Clippy retains the
existing `uv/profile.rs` argument-count warning.

The projected-4D retry reaches a 66,255,441-byte lowered scalar, matching the
baseline size and count of 47,202 family calls. Retained-body collection leaves
66,190,997 bytes including the existing definition. Sharing then creates
272 aliases in 0.245 s and reduces stored size to 7,972,529 bytes. Evaluator
expression preparation takes 0.220 s and produces a 9,560,370-byte builder input,
compared with 66,255,467 bytes previously. Thus the missed body is now compacted,
but the compact input still does not yield a completed evaluator.

Assembly completes at 536.645 s; evaluator preprocessing takes 74.671 s; Symbolica
building starts at 612.128 s. The process again crosses the 12 GiB RSS cap and is
stopped at 679.155 s, with sampled peak 12.455 GiB, after about 67 s in the builder.
The previous run stopped at 790.598 s after about 72 s in the builder. The earlier
finish of this capped attempt is largely upstream assembly variation, before
this patch acts; it is not a successful evaluator timing or a demonstrated
end-to-end speedup. Host load is recorded in the resource samples. No task build
ran concurrently; a three-second profile and subsequent report processing did.

The builder profile contains 1,186 userspace cycle samples over 2.836 s, with zero
lost samples. Hashing accounts for 62.57% of self cycles, `linearize_impl` for
5.48%, and lookup in an `AtomView`-to-`Slot` map for 3.73%. This identifies heavy
hashing and expression-lowering work during construction. The hot leaf samples
have no recovered caller frames. The sample does not measure retained allocations
or distinguish parameter-binding lookup from the subexpression cache, so the
exact memory-growth mechanism remains unresolved.
There is no completed state for a full-graph numerical or UV check. The benchmark
and profiler have stopped. Identities, receipts and checkpoint measurements are
in `soft-ct-compact-energies-2026-09-24/alias-body-compaction-2026-09-27.json`;
raw artifacts remain in `/tmp/soft-ct-alias-body-compaction-2026-09-27/`.

#strong[September 27 retained numerator evaluator bodies:] The previous
9,560,370-byte measurement covers the root and scalar aliases, but excludes
registered numerator function bodies. At that boundary, generated numerator
components still used Symbolica's default `Always` inlining. Different argument
tuples therefore lowered copies of their shared bodies. This is a concrete
integration choice; the earlier CPU profile alone does not establish an upstream
memory bug.

Change `spszwlyz` registers generated scalar numerator functions with
`InliningPolicy::Never`. Hyperdual evaluators retain `Always` because the pinned
Symbolica vectorizer rejects retained sub-evaluators. Taylor algebra, factor
collection, scalar sharing, optimizer settings and the Symbolica revision are
unchanged. Function-map metadata preserves policy and caller-scope aliases across
rebuilds and archives. The September 29 run used saved-state version 8, amplitude archive version 8,
and cross-section archive version 12. The rebased implementation uses versions
11, 11, and 14 respectively; older states must be
regenerated. Emitted standalone scripts also pin the matching numeric libraries.

Formatting and Clippy pass, with the existing UV profile argument-count warning.
All 96 focused tests pass. New checks cover smaller stored programs, independent
numerical values and first derivatives, eager and SymJIT execution, f64/f128/
arbitrary precision, binary rebuilds, and JSON standalone function metadata.
A standalone backend probe checks nested retained calls, captured globals and
inactive singular branches. A synthetic 128-call, 128-term comparison retains
514 instructions including the body instead of 33,025 inlined instructions,
with equal numerical values. That synthetic reduction is not a full-graph timing.

This experiment targets diagram 2 from the September 24 drawings: the
three-loop bubble-dressed double triangle. It starts from the rebased soft-CT
stack before the functional scalar-alias experiments. Automatic scalar-network
alias creation is disabled. Existing numerator coefficient families, which
preserve tensor factorization through Taylor projection, remain in place.

Thirty of the 33 diagram/route cases generated a usable evaluator. The 27
fully subtracted smaller cases passed their summed UV profiles. The three
child-only controls generated successfully and passed the matched soft-ray
check; their all-limits UV sweep intentionally includes unsubtracted sectors.
Diagram 2 timed out at one hour in all three routes, during evaluator
preprocessing. Its peak RSS was 9.914 GiB in direct 3D and 2.348 GiB in projected
4D. Compact energies preserve the required differentiation, but this experiment
does not make the full diagram practical yet.

== Representation and ordered operations

In the original September 24 implementation, an active on-shell energy used
`Esurface(owner, m^2, kx, ky, kz)`, the former name and component signature of
today's `OnShellEnergy(owner, s)`. The arguments contain the mass-squared
and three routed spatial components, including spatial momentum shifts and every
Taylor parameter. Neither a square root nor a quadratic momentum norm is stored
inside the energy call. There is no additive energy shift or per-expression
definition table. The owner is retained for enclosing UV mass rearrangements.
Terminal one-argument `OSE(edge)` parameters retain their existing meaning.

The semantic derivatives are `dE/d(m^2) = 1/(2E)` and `dE/dki = ki/E`;
differentiation of a routed momentum or Taylor parameter follows the chain rule
through these explicit arguments. The owner is metadata and has zero derivative.
Numeric signs in each scalar spatial argument are canonicalized using energy
evenness; this distributes numeric coefficients only, not graph numerator
products. The definition `sqrt(m^2 + kx^2 + ky^2 + kz^2)` enters the evaluator's
function map at the numerical boundary.

+ Build the routed direct-3D energy functions and collect complete factors.
+ Apply the existing hard/soft deformation, keeping active component ownership
  and both UV mass roles. Normalize mass-squared by the square of the common
  rescaling parameter and each spatial component by that parameter, then extract
  one power of it from the energy.
+ Include the branch's existing spatial measure factor `r^(3L)`, then replace
  `r` by its inverse before expanding at zero. Numerator definitions carry no
  extra measure factor. The weighted Laurent endpoints remain zero for ordinary
  U and minus one for soft S.
+ Compute the required Taylor/Laurent coefficients using the custom derivative.
  The existing numerator-family boundaries protect graph numerator products.
+ Recover Horner forms and common factors after series construction and after
  removing the auxiliary rescaling parameter, then again after forest assembly.
+ Expand the existing spatial dot-product notation at the final normalization
  boundary and perform finite tensor-component contractions. Collect scalar
  factors before evaluator construction.

The custom derivative still uses Symbolica's native function-series algorithm.
That path differentiates successively and expands a coefficient after substituting
the expansion point. It can therefore construct large temporary scalar sums
before the surrounding Horner/factor collector regains control. The compact
energy function does not by itself bound those internal derivatives. A stack
sample in the first component-energy diagram-2 run found this derivative path
growing an additive atom; sampled process RSS reached roughly 7 GiB during
that portion of generation. The pre-series checkpoints are batched: a subsequent
stage prepared 105 atoms totaling 5,250,764 bytes, with a largest atom of 202,231
bytes. The last logged atom is not necessarily the active derivative input.
The sampling pause is recorded separately with the timing artifacts.

A later evaluator-preparation sample reached `combine_numerator_families`, then
Spenso's `collect_with_map`, Symbolica's `collect_symbol`, and
`to_polynomial_in_vars`. It was cloning the polynomial variable map, including
an atom with a 101,938-byte payload. This sample precedes tensor contraction and
is distinct from the earlier derivative bottleneck. Existing temporary
coefficient aliases inside Spenso's collector are resolved before returning;
this experiment disables retained automatic scalar-network aliases, not that
pre-existing internal implementation detail.

The collector does not distribute graph numerator products. It temporarily
protects whole functions and powers, plus tensor sums inside products. Inverse
denominators and branch guards retain occurrence-local protection, so a singular
denominator cannot move outside its selector. These wrappers embed their input;
they are not retained aliases. Each accepted Horner/factor-collection iteration
strictly decreases packed atom size. Energy powers are not replaced by their
quadratic arguments during nested Taylor operations: that would erase ownership
needed by an outer U operation. Only the separate final evaluator copy removes
that distinction: equal energy arguments use a common owner, and even integer
energy powers reduce to powers of the quadratic argument. The Taylor-side atoms
retain their original ownership.

== Diagram correspondence

#table(
  columns: (auto, 1fr),
  [Drawing], [Timing fixture],
  [1], [Two-loop double triangle, ordinary-U control],
  [2], [Three-loop bubble-dressed double triangle, full soft-CT forest],
  [3], [Appendix B.1 nested gluon self-energy],
  [4], [Massive-top bubble vertex, child-only and full forests],
  [5], [Nested top self-energy],
  [6], [DGSE],
)

The remaining fixtures are the one-loop quark and gluon self-energies, and
scalar graphs with disconnected or consecutive soft components.

== Timing protocol

The numerical routes are localized 3D, explicit 3D, and local 4D projected to
3D. Raw 4D parametric integrands are explicitly unsupported by the orchestrator
and are not reported as numerical evaluator timings.

Each case starts with a fresh state, all orientations, disabled integrated and
threshold counterterms, one graph-generation worker and eight Rayon threads.
Generation and UV profiling are separate processes. The profile uses seed 1337,
25 points, scaling exponents 8 through 12, and all selected limits; the localized
route also records individual orientations for the smaller fixtures. Diagram 2
checks the full sum of 102 orientations first. Wall time, CPU time and peak
process-tree RSS are measured. A sampled 12 GiB RSS guard bounds every run. The initial
diagram-2 attempt used a 900-second stage limit; subsequent diagram-2 runs allow
3600 seconds, while smaller fixtures retain 900 seconds. A stopped run gives a
lower bound rather than a completed timing. These are single runs on a shared
host, not a controlled throughput comparison.

The suite includes all six drawn topologies, both child-only and full top-bubble
forests, and the one-loop and scalar soft-projection fixtures. Generation
completion and UV-fit acceptance are reported independently. Accepted nonvanishing
UV fits require slope at most -0.9 and R-squared at least 0.99; vanishing missing
fits retain the existing CLI acceptance policy.


== Correctness checks and existing limitations

The selected approximation, numerator and evaluator test set ran 199 tests:
194 passed and five failed. All tests specific to compact-energy differentiation,
final energy identities, direct-3D soft projection and evaluator scalar handling
passed. The five remaining failures also reproduce in the earlier soft-CT test
binary, whose checksum is recorded with the verification artifacts:

- Three paper/tensor comparisons retain a nonzero symbolic color difference:
  `paper_appendix_b_10_b_12_massless_child_and_local_b_16_nested_branch`,
  `paper_appendix_b_6_b_7_production_outer_h2_matches_the_routed_rhs`, and
  `production_outer_u_grades_a_free_fermion_mass_but_keeps_muv_fixed`.
- `integrated_independent_mixed_mass_tensors_keep_laurent_cross_terms` fails in
  FORM on a printed complex coefficient. Integrated counterterms are disabled
  in this timing campaign.
- `direct_nested_scalar_banana_replays_complete_cff_at_depth_three` has the same
  approximately 24.6% direct/projected discrepancy at depth three, stage two,
  in both binaries. This experiment does not resolve that existing discrepancy.

The initial two-argument, owner-retaining diagram-2 localized attempt completed symbolic
orchestration in 654.7 seconds, then hit its 900-second total limit during
numerator/evaluator preparation. Its measured peak was 3.47 GiB. It did not
produce a completed evaluator or numerical UV verdict. That attempt is retained
separately from the final component-argument timing matrix. A second
two-argument attempt was stopped after a stack sample found the existing final
`expand_dots` / `normalize_dots` path working through the large quadratic
arguments. The final five-argument representation avoids storing those norms
inside each energy function.


== Component-energy diagram-2 measurements

Both direct-3D routes reached the 3600-second wall limit during evaluator
preparation, before a completed evaluator. Localized 3D took 3600.684 seconds
with peak RSS 9.914 GiB; explicit 3D took 3600.648 seconds with peak RSS
9.914 GiB. Symbolic integrand construction completed in 1563.273 seconds and
1568.258 seconds respectively. No numerical UV profile was produced for either
attempt. These are censored timings, not evidence of successful full-diagram
subtraction. The localized wall time includes 25.205 seconds spent in the two
recorded debugger attachments; explicit 3D had no debugger attachment.

The projected local-4D route completed symbolic integrand construction in
529.940 seconds and entered evaluator preprocessing. It also reached the
3600-second wall limit, at 3600.634 seconds, with peak RSS 2.348 GiB. It did not
produce a completed evaluator or a numerical UV profile. No debugger was
attached to this route. All three attempts stopped with the last stage
`evaluator_stack_parse_atoms_start`; only the localized route was sampled to
identify the polynomial variable-map copying inside numerator-family collection.


== Child-only control and soft-ray check

The child-only top-bubble card sets other UV sectors to `Unsubtracted`; it
selects only the massive-top child on edges 7 and 8. The uniform all-limits UV
sweep resolves 13 of 27 summed limits as passing and labels 14 as unstable fits,
with fitted slopes close to zero, identically in the three routes. Localized
3D has 400 of 810 orientation/limit pairs passing. These failures are retained
in the timing matrix; this card does not represent full UV subtraction.

The separate matched soft check follows the existing integration-test criteria:
`top_bubble_vertex S(e6)`, seed 1337, 25 points, exponents -2 to -5 and canonical
LMB edges 6 and 8. Ordinary U and H have identical routed-ray fingerprints.
All three summed routes give U slope -2.001104449 and H slope -0.000573379,
with R-squared above 0.99999997: an improvement of 2.000531070 powers. All 30
localized orientations satisfy the soft criteria, including the eight divergent
ordinary-U orientations. Soft-check wall times and exact commands are recorded
separately in `child-soft-checks.json`; they do not alter the 33-case matrix.

#pagebreak()
== Full timing matrix

All 33 bounded cases finished on September 25. Thirty generated successfully;
27 passed the full summed UV check, and the three child-only controls have the
scope limitation above. Every fully subtracted localized smaller fixture also
passed all per-orientation UV checks. All completed generation runs recorded
zero automatic scalar aliases. The diagram-2 runs stopped before that logging
checkpoint, so their measured alias counts are unavailable.

// Summary measurements only; the raw run receipts are not source assets.
#set text(size: 9pt)
#table(
  columns: (1.05fr, 1fr, 1fr, 1fr),
  inset: 3pt,
  align: left,
  table.header([Fixture], [Localized 3D], [Explicit 3D], [4D projected to 3D]),
  [2: Double triangle + bubble], [G: 3600.68 s stopped\
Peak: 9.914 GiB\
P: —\
UV: not completed], [G: 3600.65 s stopped\
Peak: 9.914 GiB\
P: —\
UV: not completed], [G: 3600.63 s stopped\
Peak: 2.348 GiB\
P: —\
UV: not completed],
  [Quark self-energy], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.47 s\
UV: 2/2 (pass)], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.42 s\
UV: 2/2 (pass)], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.47 s\
UV: 2/2 (pass)],
  [Gluon self-energy], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.5 s\
UV: 2/2 (pass)], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.47 s\
UV: 2/2 (pass)], [G: 0.62 s\
Peak: 0.042 GiB\
P: 0.47 s\
UV: 2/2 (pass)],
  [6: DGSE], [G: 3.62 s\
Peak: 0.082 GiB\
P: 6.32 s\
UV: 27/27 (pass)], [G: 3.22 s\
Peak: 0.077 GiB\
P: 0.72 s\
UV: 27/27 (pass)], [G: 0.92 s\
Peak: 0.075 GiB\
P: 0.51 s\
UV: 27/27 (pass)],
  [5: Nested top], [G: 2.62 s\
Peak: 0.087 GiB\
P: 3.47 s\
UV: 21/21 (pass)], [G: 2.57 s\
Peak: 0.065 GiB\
P: 0.57 s\
UV: 21/21 (pass)], [G: 0.92 s\
Peak: 0.051 GiB\
P: 0.57 s\
UV: 21/21 (pass)],
  [Disconnected scalar], [G: 3.72 s\
Peak: 0.063 GiB\
P: 0.97 s\
UV: 12/12 (pass)], [G: 3.57 s\
Peak: 0.056 GiB\
P: 0.62 s\
UV: 12/12 (pass)], [G: 0.82 s\
Peak: 0.044 GiB\
P: 0.57 s\
UV: 12/12 (pass)],
  [Consecutive scalar], [G: 3.77 s\
Peak: 0.06 GiB\
P: 0.92 s\
UV: 12/12 (pass)], [G: 3.57 s\
Peak: 0.058 GiB\
P: 0.62 s\
UV: 12/12 (pass)], [G: 0.77 s\
Peak: 0.044 GiB\
P: 0.62 s\
UV: 12/12 (pass)],
  [4: Top bubble, child], [G: 3.47 s\
Peak: 0.071 GiB\
P: 52.88 s\
UV: 13/27 (failed)], [G: 3.22 s\
Peak: 0.071 GiB\
P: 1.62 s\
UV: 13/27 (failed)], [G: 1.27 s\
Peak: 0.112 GiB\
P: 2.27 s\
UV: 13/27 (failed)],
  [4: Top bubble, full], [G: 9.48 s\
Peak: 0.106 GiB\
P: 61.49 s\
UV: 27/27 (pass)], [G: 8.78 s\
Peak: 0.095 GiB\
P: 1.5 s\
UV: 27/27 (pass)], [G: 1.97 s\
Peak: 0.171 GiB\
P: 2.37 s\
UV: 27/27 (pass)],
  [3: Appendix B.1], [G: 9.88 s\
Peak: 0.348 GiB\
P: 13.43 s\
UV: 21/21 (pass)], [G: 8.98 s\
Peak: 0.349 GiB\
P: 1.17 s\
UV: 21/21 (pass)], [G: 5.42 s\
Peak: 0.442 GiB\
P: 2.04 s\
UV: 21/21 (pass)],
  [1: Double triangle control], [G: 7.32 s\
Peak: 0.098 GiB\
P: 8.58 s\
UV: 24/24 (pass)], [G: 5.97 s\
Peak: 0.084 GiB\
P: 0.77 s\
UV: 24/24 (pass)], [G: 2.22 s\
Peak: 0.065 GiB\
P: 1 s\
UV: 24/24 (pass)],
)

#set text(size: 10pt)

G is generation wall time; P is UV-profile wall time, both including process
startup. Peak is generation RSS. “Stopped” is a limit, not a completed timing.
UV counts accepted summed limits; CPU time and stage timings are retained in
`timings.csv`. Detailed orientation-level JSON results are not distributed.
A profile may contain failed fits even when its process exits successfully.


#pagebreak()
== Diagram-2 memory and reproducibility

#figure(
  image("soft-ct-compact-energies-2026-09-24/diagram2-memory.png", width: 100%),
  caption: [All three diagram-2 generation attempts. Each stopped at one hour.
  The direct-3D peaks occur during symbolic construction; all routes then spend
  substantial time in evaluator preprocessing.],
)

The generation binary is the dev-optim build with `ufo_support`, Symbolica
revision `6cdca414400ac6201ff6f692a5fa9ff21da5da58`, and SHA-256
`41af23019c808eebf1c674ee9856ce0fb2696ae799eb4bf94651f6675697bbb9`.
The jj change is `yllzpryx`; its build-time commit is
`ddfc8d858a1c3a3acc8d8970b367ef29e5151857`. Source hashes and benchmark runner
hashes were checked against the build manifest after the campaign. The final
change also contains the measurement artifacts and this report.

The adjacent `soft-ct-compact-energies-2026-09-24/` directory retains
`timings.csv`, the memory figure and run cards. Detailed JSON results, memory
traces, command plans, raw states and logs are not distributed. The original
campaign directory was `/tmp/soft-ct-compact-energies-2026-09-24/`.
Timings were collected sequentially on shared host `itphlies` with eight Rayon
threads and one generation worker. They are single observations, not a
controlled speedup measurement against main or the alias implementation.

The measurements isolate two remaining costs: large temporary scalar
coefficients inside the native function-series algorithm, and numerator-family
preprocessing with polynomial-variable-map copying. Horner/factor collection
at surrounding boundaries does not control either internal operation. An
energy-specific series recurrence and a collector that shares its variable map
are possible follow-up experiments; neither is part of this change.

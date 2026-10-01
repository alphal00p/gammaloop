#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)
#set heading(numbering: "1.")
#set par(justify: true)

= Investigation of all 34 soft-CT NixCI failures

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

Tested source: `652910c37a12981098c2f41a373de70eda93b5d5`.
Symbolica: `578dfcb55fb0662d871456ea8ecc57c1390caa77`.
Investigation: October 1, 2026, jj change `kxzkkqsm`.

The complete captured NixCI results contain 2,942 passing tests, 34 failures,
and 308 skipped tests. The sections below summarize the failure categories, differing boundaries,
evidence, limitations, and recommended corrections. The original per-test
JSON inventory and logs are not distributed.
Thirty-three identities appeared in the earlier baseline; matching names alone
do not establish matching causes. The integrated-spectacles cancellation panic
is the additional failure.

This change records an investigation. Production code, test assertions,
snapshot expectations, dependency pins, and documentation attestations are
unchanged. Diagnostic continuation in disposable fixtures is labelled explicitly;
it does not turn the original red CI results into passes. No graph numerator
was expanded for diagnostic normalization.

#table(
  columns: (auto, 1fr),
  table.header([Failures], [Established cause]),
  [15], [Idenso snapshot presentation drift],
  [7], [Documentation prose policy and stale reviewed code scopes],
  [7], [Test representation and exact-normalization mismatches],
  [3], [Malformed complex coefficients in the FORM adapter],
  [1], [Floating coefficients admitted to rational cancellation],
  [1], [Non-conserving subgraph momentum routing],
)

== Production failures

=== Rational cancellation accepts floating coefficients: one test

`scalar_spectacles_integrated_matches_across_local_uv_routes` fails while
building the first localized-3D route for the three-loop overall-log spectacles
graph. It does not reach the numerical route comparison. The parent Taylor
series completes; the subsequent optional scalar-pole cancellation panics.

The first floating coefficient appears across the Vakint evaluation boundary
for a one-loop squared vacuum propagator with exact numerator 1. The retained
child includes an epsilon-squared term numerically equal to
$-i zeta(3)/(48 pi^2)$. The decimal matches the registered AlphaLoop constant at
100-digit precision. Its origin at Vakint input/output is directly traced;
identification of the internal AlphaLoop substitution additionally uses source
inspection and the numeric match, since no internal backend-selection event
was emitted.

The one-, two-, and three-loop controls request two, three, and four Laurent
terms respectively. The third case retains the child through epsilon squared,
where the floating zeta constant enters. These positive epsilon orders are
required: a later pole can multiply them into a finite contribution.

`ScalarRationalSize::estimate` accepts every numeric leaf. Its size budget
therefore admits Float coefficients before `scalar.together().cancel()` attempts
rational-polynomial conversion. An isolated pinned-Symbolica probe using
`c*(x^2-1)/(x-1)` gives:

#table(
  columns: (1fr, auto, auto),
  table.header([Coefficient], [`together`], [`cancel`]),
  [Exact 1/4], [Succeeds], [Succeeds],
  [Exact i/4], [Succeeds], [Succeeds],
  [Float 0.25], [Succeeds], [Panics],
  [i × Float 0.25], [Succeeds], [Panics],
)

The correction belongs in the optional cancellation helper's domain eligibility
or fallible conversion. The check must traverse arithmetic, including denominator
and integer-power bases, while keeping function atoms opaque. Preserve unsupported coefficients; do not rationalize
approximate transcendental constants or discard required Laurent orders.
The helper was added in commit `c1d24f8e0dd19f3f0f714d9894f34b843fac5ff2`.
The recorded failure is a caller-domain regression, not a failure of Taylor
series or evidence that imaginary rational arithmetic is broken.

=== Vakint emits invalid FORM syntax: three tests

The affected tests are `scalar_pole_part`,
`integrated_independent_mixed_mass_tensors_keep_laurent_cross_terms`, and
`local_ir_integrated_generation_succeeds_across_local_uv_routes`.
They fail before their final physics comparisons.

Vakint's numeric-leaf serializer uses canonical Symbolica text and then replaces
the imaginary glyph with `*i_`. The printer omits an explicit unit numerator:
`i/16` consequently becomes `(*i_/16)`, which FORM rejects. The same issue affects
`(*i_)`, `(-*i_)`, and `(1+*i_)`. Actual FORM 5.0.0 rejects all five unit-imaginary
cases in the eight-case probe; nonunit controls `2i`, `3i/16`, and `1+2i` succeed.

The scalar-renormalization and nested local-IR cases expose the route without
projected tensor-integral preprocessing. The mixed-mass child contains unit
imaginary fractions. None establishes a Laurent cross-term disagreement:
backend parsing fails first. The printer is byte-identical at the previous
Symbolica pin, so this is an existing backend-adapter incompatibility, not a
new algebra change in the latest Symbolica update.

Correct the complex coefficient serializer at the Vakint-to-FORM boundary and
validate its emitted expressions in FORM, then rerun all affected routes.
Later assertions remain unverified behind this failure.

=== Nested banana routing: one test

`direct_nested_scalar_banana_replays_complete_cff_at_depth_three` shares the
physical scalar graph, unit local numerators, original-graph parameter builder,
mass substitutions, and numerical evaluator. Direct full-CFF Taylor and
projection of the 4D sectors choose different enclosing subgraph bases.

#table(
  columns: (auto, auto, 1fr),
  table.header([Depth], [Stage], [Relative route difference]),
  [2], [1], [0],
  [2], [2], [Approximately 10^-301],
  [3], [1], [0],
  [3], [2], [0.24595056],
  [3], [3], [0.24495629],
)

At depth three, stage two, the direct terminal energy uses
$|q_1+q_2|^2$, while projection uses $|q_3-p|^2$.
Physical conservation is $q_1+q_2+q_3+q_4=p$. Thus the second radicand minus
the first is exactly
$2(q_1+q_2) dot q_4 + |q_4|^2$,
which is generically nonzero. Identical mass terms cancel in this witness.
The depth-two control has no surviving spectator $q_4$ and the analogous
identity holds. This is a symbolic denominator-routing discrepancy, not
floating-point noise or an evaluator discrepancy.

The captured stage-two sub-LMB lists physical edge 4 as external but its
routing row says `Q2 = -Q1-Q3+Q5`. Embedding that row in the full physical
chart leaves residual `Q4`. `Q5` is a literal physical edge momentum, not a
combined boundary variable. Both routes can inherit this upstream problem;
neither complete counterterm is independently certified correct here.
The independent certificate and source ownership are retained in
`evidence.json`.

The exact integer-row certificate locates the defect before any hard/soft
split: the paired-crown branch of `lmb_impl` (graph/lmb.rs:997–1011) receives
both endpoints but discards the sink and emits only source-to-root external
flow. Stage one already loses `Q3+Q4`, hidden by its order-zero result.
At stage two the wrong row violates the two vertex incidence equations by
opposite `Q4` residuals. Any split of that row gives the same nonzero error
at Taylor parameter one.

The repair belongs in sub-LMB external-flow construction. Certify all routed
rows including both ends of paired crown edges, then compare the hard/soft
split. Existing downstream reconstruction certificates validate the supplied
frame; they do not validate its original physical embedding. No conclusion
about a CFF, contour, residue-sign, or Symbolica bug follows from this test.

== Test representation failures

Two export tests redundantly JSON-decode attributes already decoded by the
Linnet DOT parser. The minimal parser probe reproduces the exact column-6 and
column-4 errors. Direct attribute use and canonical DOT reserialization pass.
In the provenance test, the production accessor's exact semantic equality
assertion already passes before the redundant test-only decode fails.
The later mass assertions in the other original test remain blocked.

The collective three-component test computes zero on the production side.
Its independent Taylor oracle is also exactly zero: rational cancellation
certifies this without expanding a graph numerator. `collect_factors` alone
does not provide a complete zero test for that rational sum.

The disconnected-soft-union test compares retained numerator calls from
separate replays with distinct occurrence scopes. The 4D product and provenance
checks already pass. Resolving both sets of stored definitions with the
existing argument-aware `Integrands::resolved()` makes every forward/reverse
3D replay and residue-map assertion pass in the disposable diagnostic.
The correction belongs at the semantic comparison boundary.

The Appendix B.6/B.7 H2 test builds its independent oracle from an unreduced
color numerator; production uses `color_simplify`. Matching that input
representation preserves the independent derivative construction and lets
the complete test pass diagnostically.

The Appendix B.10/B.12/B.16 and free-fermion-mass tests have the same initial
color-preparation omission. Aligning it exposes later representation differences.
For B.16, replacing the contracted momentum–gamma pair by its defining slash
makes the coefficients of both independent slashes and the constant separately
exactly zero. The common graph numerator remains an opaque spectator.
For the mass test, the residual is a bispinor metric acting as identity on a
gamma matrix; a standalone Idenso contraction proves that local identity exactly.

These probes establish the reached mismatches as normalization gaps, without
weakening mass grading or the paper identities. The full composed comparison
helper did not expose these identities automatically. The two diagnostic test
runs still stop at those assertions, so their later assertions remain untested;
only the B.6/B.7 color-aligned paper test was shown to pass in full.

== Idenso snapshots: fifteen tests

All fifteen failures are presentation changes. Seventeen changed payload rows
represent identical pinned-Symbolica atoms without expansion. Dropping redundant
`1` before the imaginary unit also changes where the coefficient sorts among
rendered product factors. Ordered trace/chain arguments, signs, tensor indices,
network topology, and inverse momentum factors are preserved.

An exact archived Nix binary was rerun with snapshot failure bypassed and
snapshot updating disabled, solely to inspect later assertions. No later
failure occurred; `parse_val`'s exact network-to-expression assertion also passed.
`infinite_execution` completed, so its name does not indicate a hang here.
This is diagnostic continuation, not a passing original test run.
Review and update the demonstrated snapshot strings, then rerun normally.
The printer file is identical at Symbolica revisions `70375b9` and `578dfcb`;
these expectations were already stale before the latest update.

== Documentation inputs: seven tests

Six tests stop at the prose inventory: 23 authored Markdown files violate the
canonical Typst policy. The provenance test uses top-level directory symlinks,
which prose enumeration does not descend into, and instead stops at a stale
verified code-scope digest. Across 49 scopes, nine entries referring to eight
files are stale. Every anchor exists; all relevant inputs match the tested
checkout and Nix runtime source exactly. Nix filtering is excluded.

The stale files cover state persistence, process evaluation, UV orchestration,
Spenso parsing/network construction, Vakint, and the Nix test-source rules.
One concrete stale claim is already identified: architecture documentation
says state format version 7, while the implementation uses version 8.
Refreshing checksums alone would incorrectly attest to that old statement.

A disposable boundary matrix reproduced 0/7 passing with original inputs,
2/7 after prose exclusion, and 7/7 after additionally bypassing stale digest
checks. This isolates prerequisite input failures; it does not validate those
attestations or authorize weakening policy. Migrate the 23 notes to Typst,
review the changed code scopes and corresponding prose, then regenerate
reviewed digests and rerun unchanged tests.

== Evidence and validation limits

This report retains the diagnoses and selected probe results. The original
per-test inventories and raw campaign logs are not distributed.

The original CI failures remain unfixed. Diagnostic representation alignment,
snapshot bypass, and disposable documentation checksum bypass are each recorded
as such. The investigation supplies the repair boundaries and distinguishes
exact identities from still-blocked downstream assertions. It does not certify
all soft-CT routes or the complete integrated amplitude.

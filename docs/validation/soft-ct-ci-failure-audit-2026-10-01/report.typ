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

At the October 1 investigation snapshot, production code, test assertions,
snapshot expectations, dependency pins, and documentation attestations were
unchanged. Later repairs and certified target corrections are recorded below.
Diagnostic continuation in disposable fixtures is labelled explicitly; it does
not turn the original red CI results into passes. No graph numerator was
expanded for diagnostic normalization.

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

At the 1 October snapshot, the original CI failures remained unfixed. Diagnostic representation alignment,
snapshot bypass, and disposable documentation checksum bypass are each recorded
as such. The investigation supplies the repair boundaries and distinguishes
exact identities from still-blocked downstream assertions. It does not certify
all soft-CT routes or the complete integrated amplitude.

== Reference and regression follow-up, 5 October

The subsequent fixes distinguish stale targets from implementation defects.
An independent energy-contour calculation certifies the scalar-box point.
The photon point also depends on the finite-precision threshold center, whose
constraint ordering changed; independent residues certify the current center.
The original box cache test and photon point assertions pass after these
target corrections. Their numerical gates remain unchanged.

The original two-loop photon point also passes in full in 15.9 seconds.
An independent energy-order sum with explicit Dirac matrices certifies its
previously missing snapshot. The complete test is restored to CI with the
same assertion.

The complete scalar sunrise and three-mass Mercedes tests also pass with their
original million-sample budgets and compatibility gates. Their references were
already correct. Existing momentum-basis channels repair their sampling tails
without changing the integrand. Basketball and extra-loop Mercedes instead had
references for different external momenta or UV masses. Independent integrated
forests reproduce the old inputs and provide targets for the current cards;
these corrected references require separate native acceptance runs.

The complete AA photon histogram test passes all four original million-sample
runs in 694.5 seconds. Its high-energy reference rows had reused the low-energy
inputs; independent box integrals provide the corrected rows. The original
budgets and compatibility gates are retained when restoring this test to CI.

The photonic benchmark passes both one-loop cases before reaching its
15-minute whole-body cap. Its two-loop case receives only 233 seconds and
does not finish its first 100,000-sample iteration. This run does not measure
the full two-loop runtime. Extra-loop Mercedes passes all eleven first-mass
UV checks, then reaches the cap at 320,000 of the required million samples.
Its interim estimate is compatible with the independently derived reference,
but neither capped body is a complete pass. The Mercedes card sends almost all
samples directly to quad precision through its existing escalation threshold.

Actual repairs include preserving the current half-edge owner during sewn-node
traversal, reducing a separated color-generator contraction, and parsing formal
adjoint dimensions without losing their defining expressions. Directed
self-loop symmetry factors and closed-ghost signs have independent Wick and
Grassmann certificates. An earlier repaired snapshot passes all 1,923 selected core tests. The excluded generator fixtures still need their complete grouped
inventories checked; successful early assertions do not certify later stages.

Five incorrectly separated three-loop DIS pairs need the open fundamental
color network reduced to its scalar coefficient. A Schur reduction and the
long-trace Fierz identity address this missing algebra. Independent exact
color and finite Dirac-component contractions certify the affected ratios and
the complete 191-class target. The singlet eligibility tests reject free adjoint slots and unknown matrix
operators. The subsequent equivariant rule admits exactly one free adjoint
slot, while retaining the unknown-operator and additional-slot guards.

All 277 Idenso tests pass on the earlier singlet repair, including its
eligibility tests. The later one-adjoint rule adds two more controls. A fresh native
DIS probe reproduces all six independently certified pair ratios. A separate
generator failure first diverges when the comparison skeleton treats stored
ghost arrows as physical flow. Reversing stored endpoints together with the
ghost PDG must leave the chart unchanged. The correction is confined to the
temporary skeleton; regression coverage also checks the second canonicalization
pass used by cross-section generation.

That native run preserves all physical-key weights and reduces 165 groups to
164, but three charge-conjugate pairs still split. Their first unequal boundary
is the selected left/right transformation, which exchanges ghost and antighost.
For neutral adjoint Faddeev-Popov fields, charge conjugation transforms the
ghost and antighost separately and preserves ghost number; see
#link("https://arxiv.org/pdf/2205.02012")[Eq. 22]. The correction retains this
ghost flow while conjugating quarks. It is restricted to the model's neutral
color-adjoint sector and includes the captured triangle pair as a regression.
A fresh capture and an independent signed inventory agree: 481 graphs before
CP reduction become 356 raw graphs, five are color zeros, and ten numerator
classes cancel exactly. The remaining 161 groups have signed total `-1287/14`.
The complete original test body passes all thirteen assertions in 8.63 seconds
and is restored to CI.

An earlier DIS run passes its first eleven assertions, including the
certified three-loop comparison, then reaches another stale ghost-sign target.
Independent open-color projections retain 225 nonzero graphs in the case with
an additional incoming gluon and give weight sum `-97/2`. Omitting the closed
ghost sign reproduces the old `-169/2` exactly. The corresponding case with four incoming particles
has 1,428 nonzero graphs with sum `590/3`; omitting the ghost sign
reproduces its old `1298/3`. Both raw targets are corrected. The four-incoming, one-loop case also
has a stale ghost-sign reference: the raw total changes from 36 to 28, and
the grouped total from 40 to 32. Independent color and Wick calculations
and complete native inventories agree on both changes, with the graph counts
unchanged at 62 and 54. The two-loop grouped variants expose a further algebra defect: native comparison misses six
reversal pairs in the case with an incoming gluon. Independent color and Dirac
identities predict 181 groups with total `-562/7`. Their color operator has one
free adjoint slot, so the existing singlet Schur reduction does not apply. The
guarded generator projection supplies this missing reduction. A fresh complete
native inventory agrees with all 200 independently derived classes, including
19 cancellations and all 25 accepted comparisons. It produces the predicted
181 groups and total `-562/7` in 14.45 seconds.

The grouped case took 470 seconds; 437 seconds lay in eager conversion of five
sampled scalar expressions to polynomials. The sampler also used incoming
wavefunction keys for outgoing states, leaving their spinor components symbolic.
Correct outgoing keys and lazy exact conversion address these two comparison
boundaries while retaining the sample count and factorized graph numerator.
The exact zero mask also preserves the prior unit-ratio convention for a shared
zero sample; all four original comparison and flow regression controls pass.
The grouped case with four incoming particles previously reached the 15-minute
cap. A fresh run takes 168.16 seconds and produces 1,103 groups. The independent
inventory predicts 1,087 groups, singleton weight `1700/3`, and total `11773/42`,
with all 354 historical grouped-member entries preserved. Exactly sixteen
quark-triangle reversal classes remain split: their scalar comparisons still
contain an unreduced local color identity. The graph count is an actual algebra
defect; its reference is not reset to the native result.

The adjacent cubic color contraction and its literal slot signs repair this
remaining split. Symbolica deliberately preserves the slot order of functions
whose arguments contain pattern-style underscore labels. Reconstructing such
an antisymmetric function therefore cannot infer a cyclic sign by equality;
the contraction now takes the parity directly from the three stored slots.
All 281 Idenso tests pass. The complete DIS run takes 387.05 seconds and passes
its first nineteen assertions. Its last assertion has the correct 1,087
groups, all 354 ratios, and total `11773/42`; only sixteen stored signs differ.
Each affected member contains one closed quark loop and one closed ghost loop,
so its Wick sign is positive. The old negative signs also make the fixture
expression disagree with its own numeric total by 32. The corrected whole
body passes all twenty original assertions in 310.49 seconds and returns
to normal CI with the same generation options.

The original ten-assertion ungrouped jet-generator body passes in 107.11
seconds and returns to CI unchanged. Five changed statistics targets have
independent Grassmann/Wick inventories; the other five retain their historical
zero-total references. The original grouped body also passes all nine
assertions in 477.01 seconds and returns to CI unchanged. Seven changed grouped
references have independent network and ratio certificates; the remaining two
retain their historical zero-total references.

The formal-color two-loop FullSM vacuum comparison retains 300 groups.
Independent symbolic color and signed tensor-network calculations cover all
513 raw graphs, including 36 color zeros and 38 cancelled classes. Every
retained member ratio and weight agrees with the fresh native export. Six
QCD/electroweak ratios retain symbolic `Nc`; setting `Nc=3` reproduces the
preceding numerical-color reference. Its target is corrected on that basis.
The historical formal-color ledger already contained `Nc` and retained
legacy `TR`. Independent transport identifies the same 300 physical classes
and all 128 members. With `TR=1/2`, every old ratio matches its current tensor
ratio. The singleton shift `+73/2` consists of `-3/2` from closed-ghost signs
and `+38` from removing extra directed self-loop end-flip factors.

A later sampled mode exposes an actual cancellation defect. All 444 nonzero
networks, 158 accepted merges, and member weights agree with independent
Clifford and color contractions in the graph-bound denominator charts. Nine
retained groups have exactly zero scalar coefficients: expansion alone misses
their rational cancellations. Symbolica's `cancel` reduces each summand
separately; `together` combines the bookkeeping weights over a common
denominator. The zero check uses the latter while preserving graph numerator
factorization. Both symbolic sampled modes show this same nine-group defect.
The fully numerical sampled mode instead matches all 288 independent
classes, including every ratio, weight and representative. Its corrected
snapshot keeps 288 classes, with singleton weight `59/12` and total `110/3`.
All 112 historical ratios agree after physical graph-ID transport. Correct
ghost signs shift the singleton by `-3/2`; removing extra directed self-loop
end-flip factors adds `49/2`. The physical classes and comparison mode are
unchanged. The aggregate shifts are respectively `-41/2` and `481/8`.
The original whole-body run passes its first eight assertions: six two-loop
comparisons, including both 247-class sampled targets, and both raw three-loop
comparisons. The ungrouped three-loop count remains 41,650, with signed weight
`-14295/8` rather than the historical `-1087/24`. An independent Wick enumeration
covers every connected graph orbit and confirms the new total. A separate
incidence calculation verifies every symmetry factor. Omitting ghost signs and
restoring the extra factor of two for each directed self-loop reproduces the
historical target exactly. Independent generic-color algebra also retains the
expected 38,696 color-nonzero graphs and gives weight `-46369/48`; the historical
conventions reproduce its old `-885/8` exactly. The corrected targets both pass
their original native assertions. All local quartic color/Lorentz branches are
checked without expanding graph numerators.

The fourteen-assertion run reaches its 900-second cap during the ninth,
three-loop grouping assertion. The preceding eight take 631.64 seconds, leaving
only 268.98 seconds for that comparison. Its result is unproved; the timeout
does not establish a mismatch or its standalone runtime. The later filtered
assertions are not reached. Separate tests give these cases their own budgets
while retaining their original inputs and comparison settings.

The unchanged ninth assertion also reaches its own 900.26-second cap, with a
1.57-GiB observed process-tree peak and no memory abort. Its 41,650 graphs are
enumerated in 2.81 seconds and canonicalized in 5.56 seconds. The unresolved
cost is after canonicalization, in sample preparation and numerator-aware
processing; these logs do not separate parsing, tensor/color reduction and
comparison. This run predates the unchanged-flow parse-reuse guard. It provides
no completed snapshot and does not justify changing the grouped target.

A profiled replay after parse reuse also reaches 900 seconds, with a 1.99-GiB
observed process-tree peak. Basic profiling emits 3.55 million per-operation
lines as well as its 175 periodic reports. It therefore does not provide an
unprofiled timing comparison. At the last completed report, network execution
and multiplication parsing each have about 79 seconds of inclusive spans;
addition parsing has 43 seconds. These timers overlap and do not cover color
canonicalization, scalar polynomial comparison or final coefficient aggregation
separately. They cannot identify the remaining bottleneck.
The excessive per-operation traces now require the existing verbose switch;
counters, timers and slow-operation diagnostics retain their basic mode.
Spatial validation now uses an immutable walk; each color term caches its
canonical form across samples. All four focused controls
pass. Their stronger malformed-input case also exposes an existing printer
panic on `Q()`; ordinary function formatting now preserves the validation
error. Archive builds and Clippy pass. The corrected unprofiled replay
also reaches its 900.23-second cap, with a 1.85-GiB observed process-tree
peak and no memory abort. Its last completed stage canonicalizes all 41,650
graphs; it provides no later stage summary or completed grouping comparison.
The 13,975-group reference remains unchanged. These runs do not isolate a
bottleneck or establish a speedup.

The filtered raw cases also preserve their expected counts, 92 and 3,150.
Independent cycle-block and deletion-connectivity calculations agree on every
selected graph and give weights `-133/3` and `-9861/16`. No directed self-loop
survives these filters. Omitting ghost signs reproduces both historical
weights exactly, so those targets are corrected. Both raw assertions pass their original native gates, including the filtered
three-loop 3,150-graph case.

The filtered two-loop grouping target changes from 69 to 67. Independent
factorized tensor and denominator-chart calculations preserve all 92 raw
graphs and prove that six massless-Z diagrams form one scalar class. The old
snapshot split them into three classes. Reinstating that historical partition
reproduces its complete old ledger; the two extra merges explain the count
change. A separate eight-assertion test retains every original two-loop command,
flag, sample count and gate. All eight pass in 40.03 seconds, including the
92-graph raw case and the 67-group comparison, and return to CI. The unresolved
three-loop grouping assertions retain their original targets in separate tests.

Before the split, the five-assertion three-loop body reaches its final grouped
comparison and fails its snapshot in 785.97 seconds. Its raw 41,650-graph,
color-filtered 38,696-graph and filtered 3,150-graph gates all pass, as does
the selected `GL0061/GL0062` Lorentz-zero check. The final comparison returns
2,022 groups rather than 1,820, with singleton weight `-8095/48` rather than
`-1537/16`. It finishes grouping in 108 seconds, with 54 color zeros, 90 Lorentz
zeros, 1,001 merges and 37 cancellations. Its changed ledger is not accepted
as a replacement reference. Independent enumeration certifies its 3,204
filtered raw graphs and all 54 color zeros. The next boundary, the 90 sampled
Lorentz zeros and subsequent numerator classes, still needs graph-identity
reconciliation. A model-only calculation selects 90 exact zero diagrams in
the same raw inventory: each contains a two-vertex fermion trace with gamma5
and fewer than four ordinary gamma matrices. The remaining factors are
spectators. This matches the native zero count without binding its numbered
identifiers or proving its nonzero comparison classes. The test does not
retain typed graph exports, so matching
numbered graph IDs does not establish matching physical graphs. The 1,820
target was introduced in `ad0f3b8f` on 17 December 2025, replacing 2,009 while
fixing sampled zero division and polynomial caching. That history alone
does not certify either current count. This is a snapshot mismatch rather
than a timeout or memory abort. The observed process-tree peak is 6.98 GiB.

The four passing assertions now run in two independent CI tests with their
original commands, sample counts and snapshots. The raw and color-filtered
inventories pass in 511.37 seconds; the filtered inventory and selected
Lorentz-zero check pass in 25.99 seconds. Their observed process-tree peaks
are 6.97 and 0.57 GiB. Both runs use one feyngen worker and no retries.
The raw test has a 15-minute budget and occupies all nextest test slots in
its integration-check process. The filtered test keeps the default four-minute
budget. Neither reservation claims an entire NixCI runner. The unresolved
filtered grouping reference remains unchanged in its separate slow test.

Kaapo's four-loop QCD reference retained two ghost diagrams whose color factors
vanish exactly. Both have the six-vertex K3,3 color network. Exchanging two
vertices on one side gives three odd permutations of the opposite `f` tensors,
so the color factor equals its negative for generic `Nc`. Reversible dimension
cooking now allows the existing zero detector to recognize this identity.
The Lorentz numerator remains a spectator. Removing weights `-1/4` and `-1/6`
changes the surviving sum from `-25/3` to `-95/12` and the count from 91 to 89.
An independent incidence calculation matches all nine native fermion/ghost
buckets. The historical 52-graph, ghost-sign-omitted `-44/3` checks remain
unchanged. No admitted graph has a directed self-loop, so end-flip weighting
cannot explain this correction.

Basketball's norm-only stability check accepted equal-magnitude values with
opposite signs after rotation. Its scalar amplitude requires signed equality.
The signed check rejects that double-precision result and falls back to a value
consistent with the independent high-precision calculation at the existing
tolerance. Double triangle has a different failure: the card retains only
virtual cuts, whose signed sum has a nonzero inverse-cubic soft term. Its
finite integral target cannot be validated for that cut selection.

A matched Double Triangle control retains the same graph, momentum basis,
eighteen orientations, and virtual terms, and enables its two real cuts.
Their signed sum cancels the inverse-cubic and inverse-quadratic soft terms;
the remaining inverse-linear term is integrable in the tested noncollinear
soft cone. This checks the missing cancellation locally, but does not certify
global convergence or the old finite integral reference.

The finite-integral fixture now requests the inclusive selection, retaining the
same graph and both virtual terms and adding the two real cuts. Its original
point agrees with independent energy residues in double and Arb precision.
A separate covariant forest calculation gives the inclusive diagram divided by
Born as `CF*alpha_s/(4*pi)*(35/3 + 2*ln(MUV^2/s))`. At the unchanged `MUV=1000`
and `s=4`, the real integral is `1.8182851767060301e-5` and its imaginary part
vanishes. Independent review checks the finite local counterterm and the
complete-forest discontinuity. The complete original million-sample test passes its one-sigma acceptance
in 112 seconds and is restored to CI. This validates the inclusive fixture,
not the historical virtual-only integral.

A second Basketball probe exposes equal, inaccurate quad-precision results
under both rotations: about $4.15 times 10^6$ instead of $3.40007 times 10^(-4)$.
Native Arb and an independent factorized calculation agree on the latter.
The card now enables the existing zero-discrepancy escalation guard, which
rejects that quad result and reaches Arb. An accurate negative control still
uses double precision. Equal accurate probes can also escalate, so this repair
can increase runtime; it does not establish convergence of the original
million-sample integral with the ordinary sampling map.

An isolated four-channel power-map trial passes the original 20 UV probes
but fails during Monte Carlo preparation after 218.68 seconds. The selected
channel has `log(J_forward*q_inverse)=0.00147714`, exceeding its unchanged
`1e-6` density budget at quad precision. Sampling fixes its source precision
before drawing proposals; later integrand precision retries reuse that proposal
and cannot repair its map. The failure occurs before statistical acceptance
and does not contradict the integral reference. The trial card is not adopted.
The error does not retain the failing coordinates, so radial inversion and
loop-basis conversion remain distinct unresolved causes.

The strict TTH point keeps its historical target. It was introduced on
15 December 2025 in `1fc2cb1b`, selecting `GL18` with CP grouping, without a
stored typed graph. The surviving TTH benchmark graph has four loops and
cannot identify that three-loop NLO fixture. Current CP-grouped `GL18` matches
current ungrouped `GL22`, but this does not establish either graph's identity
with the original. The following day's vertex-momentum permutation is excluded:
all admitted nonzero vertices at this QCD order are momentum independent.
Potential momentum-dependent gluon or ghost vertices instead produce disconnected
components or vetoed zero tadpoles. This excludes that permutation hunk, not
other changes in the same commit. Later graph-library and external-assignment
changes remain comparison boundaries, not demonstrated causes. The missing
certificate is the original typed incidence, external labels, vertex ports and
signed loop basis. Its strict mismatch remains unresolved, and changing the
target would hide that uncertainty.

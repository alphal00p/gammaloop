= Vacuum grouping B2 diagnosis
<generation-integration-wp1>

== Verdict and scope

#strong[Verdict: permitted grouping; the reported sample-dependent coefficient
is a presentation artefact.] Frozen v4 Group(28) is the numerator of GL028
relative to GL034. At the production SU(3) convention, their original routed
numerators satisfy the exact factor-preserving identity
$N_28 = (7/2) N_34$. Exact cancellation of each observed stored scalar quotient
also gives 7/2: no sample-coordinate or MT dependence survives.
This is a bounded certificate for that pair, not certification of all generation
snapshots or of the integration. Independent Astra review is required.

The implementation workspace is
`/shared/lcnbr/dev/gammaloop/generation-b2-diagnosis`, based on frozen
`e0a779808174c64057fecfd062c948b8de8e5c1f`. The upstream control workspace is
`/shared/lcnbr/dev/gammaloop/generation-b2-upstream`, based explicitly on
`588148851e3e989ecbb7478279db7d46727e23a8`; its durable bookmark is
`generation-b2-upstream-control`. No frozen reference or donor checkout was
modified. No production correctness change, WP3 port, push or PR belongs here.

== Complete pipeline comparison

Both routes request the Standard Model with its default restriction, four total
loops, no external legs, g/ghg/t/u/d, QED order zero and scalar grouping.
Their ordered stages are model restriction, enumeration/admission, physical
routing and local-factor construction, color simplification, common comparison
routing, optional canonical comparison, component sampling, scalar quotient
comparison, collapse and weight printing. Upstream combines these phases;
v4 separates generation, build and grouping.

The last shared boundary is the request, not identical actual model data.
`d8146dc2ad8c` changes UUV1 from P(3,2) to P(3,2)+P(3,3), reduced during
generation as −P(3,1). Enumeration/admission and canonical comparison frames
also change later. The full weight sums therefore do not have identical inputs;
quotient normalization alone cannot explain their difference. The selected
GL028/GL034 pair contains gluons and massive top quarks, with no ghost vertex,
so UUV1 is excluded for its own identity check.

No UV forest, CFF construction, residue, subtraction or integration stage is
reached by `--only-diagrams`. Those mechanisms are excluded. A Group index is
a source index, not a revision-independent physical identifier: upstream's
Group(28) belongs to a different labelled pair. The archived source expectations
are distinguished from the newly observed outputs.

#table(
  columns: 4,
  [*Route / mode*], [*Samples / seed*], [*Groups*], [*Cancelled weight sum*],
  [Upstream production], [5 / 3], [176], [1157/168],
  [v4 production], [5 / 3], [176], [1607/168],
  [v4 production, varied coordinates], [5 / 7 and 11], [176 each], [1607/168 each],
)

These sums are generation bookkeeping diagnostics, not amplitude-equivalence
proofs. The change of full sums is not investigated beyond its differing input
boundary in WP1. WP9 owns the remaining snapshot certification.

== Original pair and exact certificate

The annotated fixtures
`crates/feynkit-generator/tests/fixtures/vacuum_b2/GL028.dot` and `GL034.dot`
were captured at the accepted grouping comparison boundary before collapse.
They retain the original vertex and edge factors, prefactors, projectors,
incidence, local indices, signed routing, loop basis and model fingerprint.
`evidence/generation-b2/pair-provenance.json` records hashes and direction:
candidate source 28, master source 34. No numerator was replaced by a fitted
ratio or a nonlocal synthetic polynomial.

Independent Astra ROUND1 requested changes at `468d10cf`. The old proof
called `AlgebraSettings::hep()` before equality: Dirac monomial collection and
default expanding contraction distributed propagator alternatives. Its later
exact zero result did not close the factor-preserving gate. That certificate is
superseded, rather than repaired by a final `expand=false` flag. The independently
confirmed stored-quotient cancellation and fixture provenance remain valid.
See `evidence/generation-b2/review-round1.json` and the retained verifier sources.

The corrected test uses a closed-loop matrix similarity, preserving the original
local products through the complete symbolic argument:

+ Assert the original prefactor and projector are both one. Load the existing
  aggregate numerator unchanged, replace only the model coupling definitions,
  and route Q into the common formal loop basis with the existing
  `normalize_momentum_routing`. Coupling replacement does not distribute
  tensor sums; routing's `expand_num` normalizes numerical multiples of linear
  momentum sums without multiplying independent numerator factors.
+ Select color heads with `coefficient_list`. For these sources there is exactly
  one color sector. The unselected Lorentz propagator sums remain immutable
  spectators. A syntactic reconstruction assertion requires the factored
  color-times-Lorentz product to equal the original complete routed expression.
+ Contract only the isolated color product with the exact SU(3) component
  library (`hep_lib_atom`) and the finite tensor-network executor. No Dirac or
  Lorentz expression is supplied to that executor. Assert the resulting color
  constants are 7/3 for GL028 and −2/3 for GL034. Scalar cancellation here acts
  only on the completed constant color result.
+ Use the production Weyl gamma components and charge-conjugation matrix to
  verify $C^T=-C$, $C^2=-I$ and $C gamma^T C=gamma$ exactly for every matrix entry and all
  four Lorentz components. These are finite constant matrix identities; no
  momentum or mass is evaluated.
+ Assert that vertices 1,2,3 and fermion edges 3,4,6 form a closed loop. Audit
  its entire matrix domain: six four-dimensional gamma factors and three
  identity factors from the mass alternatives, with no unsupported matrix
  head. Every loop spinor slot occurs in its vertex gamma and both alternatives
  of exactly one propagator; this checks that the similarity domain is closed. Apply $T(A)=C^(-1) A C$ to each local matrix: identities stay fixed and
  each gamma becomes its negative transpose. The same similarity on every
  factor preserves the closed trace because adjacent C and inverse-C cancel.
  Linearity applies inside each intact mass-plus-slash factor; the test never
  distributes those factors or substitutes values for their scalar coefficients.
+ Normalize existing numeric content and canonically relabel dummy indices
  directly on the Atom graph. No tensor finisher, symbolic contraction,
  Dirac simplifier, `AlgebraSettings`, or algebra prerequisite is called in
  this step. Literal equality of the transformed candidate and negative master
  certifies $L_28=-L_34$, retaining all formal loop-momentum, MT and G dependence.
  Thus $N_28/N_34=-(7/3)/(-2/3)=7/2$ exactly.

The separation, component constants, matrix-domain guards, and final literal
identity are executable assertions. The component library is used only for
constant color and gamma/C tensors, not for a finite-point test of a symbolic
graph numerator. Unchanged local vertex/edge factors and signed routing keep
the locality independently verified by Astra. The proof uses no whole-numerator
polynomial expansion and no on-shell or sampled-kinematic assumptions.

The actual production control explicitly sets
`check_canonical_numerator=true` and `number_of_numerical_samples=0`.
It must return two singleton groups, no removed source, and each source 34 and
28 exactly once, with unit ratio, original name, intact numerator, valid owner
index and unchanged scalar factor/weight. Five samples must instead return one
group containing both actual source IDs, ratios 1 and 7/2, valid owner indices
and total scalar weight 3/4. The zero-sample two-group result and canonical-plus-
sampled one-group result are an asserted fallback witness for this original
pair, not a toy comparison. These production preparation calls are coverage
checks after the separate factor-preserving proof; they do not supply it.

The scalar contract permits model parameters, couplings, functions, form factors
and scalar color invariants. MT dependence alone would not violate it. Exact
cancellation of the original stored quotients, including both original slow-test
rows, has already been independently confirmed by Astra. That display result
is distinct from the corrected source-pair certificate, which awaits a new
independent review.

== Production path and regression

Production defaults are five samples, seed 3, canonical comparison disabled,
and model scalars retained. The full-request runs use those defaults and repeat
with seeds 7 and 11. The scalar quotient's spelling changes with the points,
while exact scalar cancellation consistently gives 7/2 for the selected pair.
The durable JSON retains the accepted quotient, not merely its evaluated value.

The unit regression
`grouping::tests::vacuum_group28_original_pair_identity` retains its explicit
ignored WP1 validation entry point. It first establishes
the original-pair identity independently. It then feeds the original local
factors through `group_diagrams`, with canonical comparison both disabled and
enabled, five samples, and seeds 3, 7 and 11. Every accepted stored ratio must
cancel to 7/2. At unused sample indices 7 and 13 it additionally checks each
original numerator against the accepted ratio times its original partner.
These checks supplement the exact certificate; they do not replace it.

From the workspace root, silently preserve or load `SYMBOLICA_LICENSE`, then run:
```sh
nix develop -c cargo test -p feynkit-generator --profile dev-optim --lib \
  grouping::tests::vacuum_group28_original_pair_identity \
  -- --exact --ignored --nocapture
```
The bounded `vacuum_b2` example regenerates the full-request outputs and source
provenance; its arguments are output directory, seed and sample count.
The exact commands, pinned shell realization, outcomes and retained log hashes
are recorded in `evidence/generation-b2/execution.json`; the correction and
prerequisite audit are in `evidence/generation-b2/fix-round1.json`.

No correctness defect was reproduced, so no bisect or first-bad-commit verdict
is warranted. The earlier pre-v2 anomaly and later quotient normalizers remain
historical provenance, not attributed causes of a physics error. There is no
production remedy to port in WP3 for this pair. WP3 must carry the original
fixtures and signed provenance, this corrected certificate, asserted zero-sample
control, sampled/fallback matrix, and original stored-quotient cancellation
evidence, then rerun the ignored regression on the rebuilt revision. Any scalar-only presentation
normalization or snapshot update belongs to a separately reviewed step; matching
counts or sums alone must not certify those snapshots.

== Validation limits and retention

The original slow test was exercised on both frozen source trees. Ordinary
assertions fail at the existing two-loop snapshot before reaching four loops.
A separate diagnostic run set `INSTA_FORCE_PASS=1` and `INSTA_UPDATE=no` to
execute all its requests, including the zero-bridge / one-factorized-subtopology
variant. Those runs are explicitly not passing snapshot certifications.
Their structured differences and ordinary assertion failures are retained.

A full 270-diagram v4 exact-only run was stopped after 605 seconds, at roughly
9.9 million KiB RSS, without a grouping result. This bounded control is recorded
as incomplete, not a bad verdict. The completed original-pair exact certificate
and focused production-owner controls are the correctness evidence for B2.
The full exact-mode cost remains outside this bounded diagnosis.

Retain both task workspaces, their named bookmarks, the small DOT fixtures,
structured JSON, and task-owned `target/wp1`, `target/wp1-review-round1`,
and `target/wp1-fix-round1` evidence/logs for independent review.
Temporary expanded exact diagnostic outputs and redundant source dumps are
removed after their receipts are preserved. Shared caches are untouched.
NixCI remains disabled. No background process is intentionally retained after
validation finishes. No source or historical cleanup was performed.

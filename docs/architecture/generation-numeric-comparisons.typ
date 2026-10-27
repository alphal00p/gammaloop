= Generation integration numeric comparison slice

== Scope and source mapping

This independent WP3 slice ports source units 08 and 09 onto approved prefix
`74cebf6c8f673f5849b26afe55fa98035d5eba1d`. It changes the existing
`ApproxEq` owner in `crates/gammalooprs/src/utils/mod.rs`.
The structured receipt in
`docs/architecture/evidence/generation-numeric-comparisons/receipt.json`
records exact source identities, implementation identity, commands and outcomes.
Both sources were initially ported in
`19e0e9c27642523e0a430cdba0259d380e65c4f7`. Source 08 additionally requires
precision correction `c72726011f70640f1e20495f4f62d969afef6743` after Astra round 1.
Source 09's production implementation remains unchanged by that correction.
Parent replay tracking is reconciled only after independent Astra approval.

- Source 08 (`2e77886cde8d72e1e063cafced9a2738811abf20`): compare the circular
  phase residual with the existing directed reference scale.
- Source 09 (`0eb460fd62c3b7fbd2a244ba13d8f38b60fbfad6`): compare the
  imaginary part against zero, using absolute tolerance, for complex-to-real
  comparisons; delegate the reverse direction to that same implementation.

== Invariants

Approximate equality stays in `ApproxEq`; standard equality contracts and
callers do not change. Norm and real-part comparisons retain their existing
relative scaling, including the scalar owner's directed reference scale.
Phases retain relative tolerance away from zero and the scalar owner's absolute
zero-reference rule. This slice does not introduce globally symmetric relative
comparison semantics or replace the phase tolerance by a fixed angular bound.
Imaginary parts use the strict absolute bound `abs(im) < tolerance`, independently
of the real-part scale, including when the real reference is zero.

All generic constants and intermediate phases preserve the active backend and
precision. Overflow and invalid values must remain rejected rather than treated
as underflow-to-zero. Existing norm construction and underflow semantics are
outside this slice; no tolerances or existing expectations are weakened.

== Astra round 1 response

The first independent review requested changes for P2: normalizing a native tiny
negative phase to `[0, 2π)` could round it to `TAU`, so subtracting `TAU` later
erased a representable phase residual. Both comparison APIs then accepted a
meaningful discrepancy. The retained reproducer and review receipt live in
`target/astra-numeric-round1.RJeXH7`; their recorded artifact hashes are preserved.

The correction computes the residual directly from signed native phases. Across
the negative-real branch cut, it subtracts each phase's nearby signed π endpoint
before taking their difference. Only the denominator uses the reference phase
normalized to `[0, 2π)`. A zero reference keeps the absolute residual. The existing
scalar `ApproxEq` owner compares the scaled residual with zero; there is no new
comparison API. Norm calculations and source 09's production implementation are
unchanged.

The P3 response adds exact-context `9/40` mixed-comparison controls. Complex real
part 5 against real reference 4 rejects, while complex real part 4 against real
reference 5 passes. Both mixed call directions and both signs retain this real
reference normalization; symmetry of a relative comparison is not imposed.

== Validation and handoff

The source regression intent is retained and extended to binary64, double-double
`QuadFloat`, 256-bit `SamplingFloat`, and 1000-bit `ArbPrec`. Controls cover both
comparison directions, the negative-real branch cut, the zero-angle seam,
relative phase scaling, resolvable phase differences, exact and near reals,
pure imaginary values versus real zero, strict imaginary tolerance boundaries,
real mismatch, invalid values and overflowing norms. The existing momentum
rotation test exercises the mixed scalar/complex caller.

The initial baseline failed all three original regressions; the first port passed
29 selected tests, but independent review then found P2. In the correction round,
all four backend-specific precision regressions fail against that faulty owner,
while the P3 normalization control already passes. After the correction, all
34 selected owner/caller tests pass, including the original 29 and five new tests.
The new cases use native-context decimal powers with backend-specific exponents
17, 35, 80 and 305. Both imaginary signs, both comparison directions and both APIs
are checked. Precision-scaled negative-real branch-cut controls reject and accept
at distinct fixed multiples of each backend's epsilon.

Clippy and formatting pass, with two existing warnings in vendored Clarabel.
Standalone Typst compilation, registry/title contracts and retained evidence
hashes also pass. Exact commands, result excerpts and finding responses are in
the receipt and results file. No full physics suite, performance benchmark, push
or publication belongs to this correctness slice. The workspace, named bookmark,
first Astra review/reproducer and correction logs are retained for a new
independent Astra review. Parent integration remains blocked on that approval.

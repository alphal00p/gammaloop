# Changelog

A curated, reverse-chronological record of `oneloop`'s development milestones. This is a distilled
milestone log, not a per-session diary. Entries group related work; exact numbers are quoted from the
validation record (see [benchmarks](06-benchmarks.md) and [../benchmarks/README.md](../benchmarks/README.md)).

For orientation, see [the overview](01-overview.md); for the algorithm, [the reduction algorithm](02-reduction.md);
for numerator handling, [numerators](03-numerators.md); for the current edge cases, [the frontier](04-frontier.md);
for the app integration, [the app](05-app.md).

---

## gg→h form factor — first full-amplitude validation (2026-08-22)
- **End-to-end on a physical process (`benchmarks/rust/ggh_formfactor.rs` + `python/ggh_formfactor.py`).**
  Closes the "target B" the summary deferred, for gluon-fusion Higgs (top loop). The top-loop
  Dirac trace is contracted with the transverse projector `P_{μν}=g_{μν}−(q₁q₂+q₂q₁)/(q₁·q₂)`
  (the plain `g_{μν}` keeps a non-transverse contamination that drifts the τ-shape ~1 %),
  giving `P_{μν}N = 4m(6−d)k² + 4m(4−2d)k·q₁ − (64m/s)(k·q₁)(k·q₂) + [2ms(2−d)+4m³(d−2)]`
  (verified to 14 digits vs explicit Dirac matrices). The reducer reduces each `dot(k,·)`
  monomial (exact — the on-shell external Gram `[[0,s/2],[s/2,0]]` is non-singular), OneLOop
  evaluates the masters, the `(4−d)=2ε` from the trace kills the `B0` pole, and the finite part
  gives `A_{1/2}(τ)`.
- **Numbers.** `A_{1/2}(τ)` matches the analytic `2/τ²[τ+(τ−1)f(τ)]` to **10⁻¹³** across six
  (m_H, m_t) points (physical point 1.3762612), with the exact `1/m_t` normalization confirmed
  by varying m_t. Assembled `|M|²` = `9.3663e-3` vs MadLoop `9.3703e-3` (**0.04 %**) at the
  Higgs-scale α_s; the ~12 % at α_s(m_Z) is purely the α_s scale, not loop content.
- **Tooling.** OneLOopBridge wired up for master evaluation (built via maturin against the local
  `libavh_olo.a`; py3.14 + gfortran runtime path). See [09-ggh-formfactor](09-ggh-formfactor.md).

## Reduced-numerator display compaction (2026-08-21)
- **App polish (`reduce_bridge.rs`).** The `reduced_num` was `Σ coeff × master` but each
  monomial carried a ~95-char graph-grouping/sign bookkeeping tag and was duplicated across
  grouped graphs. Now: (1) `evaluate_overall_factor` collapses `NumeratorDependentGrouping`,
  `AutG`, `InternalFermionLoopSign`, … to numbers (grouped-graph duplicates sum); (2) the
  reduction is folded into one atom and `collect_factors`'d, pulling the common
  coupling/colour/polarization prefactor out front and collecting like terms. Masters are
  emitted as opaque `A0/B0/C0/D0` symbols (flat args) so the factoring leaves them intact.
- **Effect:** `gg→h` `reduced_num` **11,875 → 1,799 chars (6.6×)**, 60 → 0 grouping tags;
  validated across A0/B0/C0/D0 (a→a, gg→h, light-by-light boxes), renders in Typst, masters
  unchanged. Worst-case 300 KB box tensor factors in ~8s (debug), well under the app timeout.

---

## Chirality-projector traces — closed fermion loops collapse (2026-08-21)
- **idenso fix (the app-numerator win).** Closed fermion loops carrying a chirality
  projector (`gg→h`, `H→γγ`, all electroweak) left an inert `Tr(γ … ℙ_p … γ)` in the
  reduced numerator. Root cause: `dirac/simplify.rs`'s trace evaluator bailed on the
  projector factor. Fix (3 surgical additions, all in `idenso`): expand
  `ℙ± → ½(δ ± γ5)` at `simplify` entry, close an equal-endpoint `chain(i,i) → trace`,
  and a `bispinor_rep_of_slot` helper. `spenso` untouched; flows to the app through the
  existing `simplify_gamma()` in `gammalooprs`; no public API change.
- **Validated end-to-end.** `Tr(ℙ₊ γμγνγργσ) →` exact metrics + ε (unit test); 21
  FORM/FeynCalc reference tests pass; real `gg→h` now reduces to
  `C0(…; MT²,MT²,MT²) + 3·B0` with **zero residual traces**. See
  [projector traces](08-projector-traces.md).

---

## Completeness + validation expansion (2026-08-21)
- **Reducer completeness.** `gram_solve` now pseudo-inverts a rank-deficient Gram (solve on a maximal
  independent sub-Gram, redundant coefficients zero), so **N≥7 tensor numerators** and coincident-momenta
  cases reduce instead of panicking — the last guarded panic in the tensor path is retired. Rank-3–6
  mixed tensors validated (box→heptagon, including dotted); the cross-engine suite is now **132/132**.
- **Bridge N≥5.** The graph→reducer bridge maps externals dynamically (`P(j)→q{j+1}`, up to nine-point),
  so the **deployed app reduces pentagons and beyond** (validated live on `e+e-→ddg`, `q4` mapped).
- **MadLoop suite → ~110 processes.** Three MG5_aMC v3.7.2 batches (96 fresh) reproduced, spanning
  triangle/box/pentagon/4-gluon/2→3-jet/loop-induced/diboson, QED+QCD, massless+massive; the IR/colour
  pole structure (massive-vs-massless double pole, exact colour factors) comes out textbook. A
  reproduction record — the direct reducer validation remains the cross-engine suite.
- **Docs.** Added the one-page [summary & document map](00-summary.md) and the
  [benchmark report](07-benchmark-report.md).

---

## gammaloop → family bridge

- Added `bridge.rs`: translates a gammaloop one-loop graph into an `IntegralFamily` the reducer can
  consume. See [the app](05-app.md) for how this feeds the FeynmanEngine integration.
- `numerator_to_dot_form(num: &Atom, heads: &GammaloopHeads) -> Atom` rewrites a gammaloop (spenso)
  contracted numerator (edge momenta carrying explicit Lorentz indices) into oneloop's `dot(k, q)` form.
- `external_offset_from_lmb_rep(lmb_rep: &Atom, heads: &GammaloopHeads) -> Atom` extracts the
  external-momentum offset `r` such that an edge momentum is `k + r`.
- `GammaloopEdge { lmb_rep, mass_sq }` and
  `family_from_gammaloop(numerator: &Atom, edges: &[GammaloopEdge], heads: &GammaloopHeads) -> IntegralFamily`
  assemble a family end-to-end: pairwise invariants `(r_i − r_j)²` computed from the offsets, propagator
  exponents set to 1, numerator translated to dot form. (`GammaloopHeads` names which symbols mean the
  loop/external momenta, Lorentz index, and metric.)
- This realizes Path A of the graph bridge (gammaloop's contracted numerator is already publicly
  emittable via dot export), so no gammalooprs API change is required for the numerator handoff.
- Tests: `bubble_invariants_from_offsets`, `scalar_massless_bubble_reduces_to_b0`,
  `rank1_bubble_numerator_reduces_through_the_bridge`.

## On-shell massless-leg regularization (`reg_delta`)

- Fixed the one remaining gap on generic on-shell processes: massless on-shell external legs (`p² = 0`)
  make pinched sub-Grams / sub-Cayleys vanish, which previously caused the reducer to panic
  (`singular Gram matrix`, `degenerate Cayley reduction: exceptional kinematics`).
- `reduce()` is now a thin wrapper. If any `kinematics.invariants` entry `is_zero()`, it swaps each
  zero invariant for a shared symbol `oneloop::reg_delta`, calls the unchanged core reducer, then
  substitutes `δ → 0` in the resulting coefficients and master arguments and drops zero-coefficient terms.
  Generic kinematics take the unchanged fast path with zero overhead.
- The fix is exact: the reducer's exact rational arithmetic cancels the `1/δ` inverse-Gram poles
  algebraically in the sum `Σ cᵢ Mᵢ`, so coefficients come out polynomial in `δ` (finite at `δ = 0`).
- Validated on the gg→h massive-top triangle, on-shell massless legs, rank-2 `(k·q₁)²`:
  reduces to `(1/20)[B0(0;1,1) − B0(2/5;1,1)] = −3.4748e-3`, matching an independent Feynman-parameter
  computation (`−3.474834e-3`) to `~3e-10`, and a numerical off-shell δ-ladder extrapolation to the
  same limit.
- Broad check across 13 configurations (2/3/all massless legs, massive vs massless internal, ranks 0–2):
  zero panics, all results physically correct, with path-independence confirmed (rank-2 vanishes as `δ²`,
  same for symmetric and asymmetric approaches). Notably the fix also handles `dot(k, q₂)` on the
  massless-internal e+e→dd̄ triangle, which previously panicked.
- Method is the exact-arithmetic analog of expanding about vanishing Gram determinants
  (Denner–Dittmaier, hep-ph/0509141) / auxiliary-mass-flow (AMFlow, arXiv:1711.09572).
- Test: `on_shell_massless_triangle_rank2_regularizes`. See [the frontier](04-frontier.md) for the
  validity domain and remaining cases.

## MadLoop validation (~20 physical processes)

- Reproduced Valentin's reference benchmark `e+ e- > a > d d~ [virt=QCD]` with MG5_aMC v3.7.2 to
  ~14 digits: Born `3.4754514769164148e-03`; virtual normalized by `Born · αₛ/2π` — Finite
  `−8.9363792407373648e+00`, single pole `+8.7724371707012754e+00`, double pole
  `−2.6666666666666670e+00` (`= −2 C_F = −8/3`), at `μ_R = M_Z = 91.188 GeV`, `s = 10⁶ GeV²`.
- Understood the pole structure analytically: double pole `−2 C_F` exact; single pole
  `C_F[−3 − 2 ln(μ²/s)] ≈ 8.77`. The loop is the massless on-shell triangle `C0(0,0,s;0,0,0)`.
- Reproduced ~20 processes spanning every one-loop topology:
  - **Triangle:** e+e→dd̄ (−8.936), Drell-Yan uū→e+e− (−8.936 by crossing), uū→Z, ud̄→W+, tt̄→g
    (massive top), e+e→bb̄ (massive b).
  - **Box:** uū→dd̄, dd̄→dd̄, uū→gg = dd̄→gg = gg→dd̄ (−54.03 by crossing), gu→gu, uū→tt̄ (massive top).
  - **Pentagon:** e+e→dd̄g (finite 2.4369).
  - **4-gluon:** gg→gg (123 loop diagrams, −66.63).
  - **Loop-induced:** gg→h, gg→hh (finite `3.4123262140874510e-05`), H→γγ [virt=QED]
    (28 W+top loops, 6.6360e-2), γγ→γγ light-by-light [virt=QED] (186 fermion-box loops, 1.0538e-3).
- Characterized the on-shell frontier definitively (see [the frontier](04-frontier.md)): triangles
  tolerate one on-shell leg; boxes are robust; scalar and rank-1 reduce on-shell; only the
  ≥2-on-shell-leg triangle at rank ≥ 2 (or a raised propagator power) hit the degenerate wall — the
  gap the `reg_delta` fix then closed.

## Cross-engine numerical validation harness

- Built the `benchmarks/` cross-check: reductions are emitted from Rust and evaluated
  against three independent oracles — OneLOopBridge (avh_olo Fortran, numeric A0/B0/C0/D0), feynalg
  (analytic masters, Denner convention), and scipy plus two Feynman-parameter tensor oracles.
- Coverage: scalar N = 3–7 (exact, generic kinematics), tensor numerators to rank 6, UV-divergent,
  massless internal lines, timelike / threshold / near-degenerate kinematics.
- Result: **110/110 benchmark families pass, 592/592 masters**, and this held after the on-shell fix
  landed. See [benchmarks](06-benchmarks.md) for the harness details.
- Fixed a robustness gap surfaced by the harness: an exactly-singular Gram (e.g. a box with `q₂ == q₃`,
  rank-deficient external momenta) previously panicked on an `.unwrap()`; it now `.expect(...)`s with a
  clear message (`singular Gram matrix: external momenta linearly dependent`). The full-rank path is
  byte-identical.

## Numerator reduction (rank-1 / rank-2)

- Added tensor-numerator reduction: numerators are polynomials in `dot(k,k)` and `dot(k, q_i)`, reduced
  via Passarino–Veltman transverse projection and FJT/Tarasov index-lowering.
- Rank-1 (linear in one `dot(k, q_i)`) and rank-2 (products of two dots) reduce for triangle and box;
  the bubble handles rank-1 dots directly.
- Numerator RSP rules interpret invariants in physics order; the dispatcher remaps
  lexicographic → physics ordering at the numerator entry point. See [numerators](03-numerators.md).

## Raised propagator powers via per-topology IBP

- Added raised-power reductions for tadpole and bubble (raised propagator exponents), folding higher
  powers back to the scalar masters through the closed-form per-topology IBP recursions.

## N-point scalar reduction (N > 4)

- Added scalar reduction for any propagator count. `reduce()` dispatches on N: 1 = tadpole (→A0),
  2 = bubble (→B0 + tadpoles), 3 = triangle, 4 = box, and **N > 4** reduces to boxes via the
  van Neerven–Vermaseren bordered-Cayley construction (with an FJT/Tarasov degenerate branch for
  `det(Y) = 0`). Validated exact for scalar N = 3–7 in the harness.

## A0 / B0 / C0 / D0 master coverage

- Established the master-integral basis: `MasterIntegral` enum with `Tadpole` / `Bubble` / `Triangle` /
  `Box` variants, all fields exact `Atom`s. The `symbol()` emitter produces opaque function atoms
  `A0`/`B0`/`C0`/`D0` in AVH/OneLOopBridge argument order. No analytic forms live in the core, by design:
  the crate is the REDUCE step; EVALUATE is delegated to OneLOopBridge (numeric) and feynalg (analytic).
- Foundation for the whole pipeline: `IntegralFamily` / `Integral` / `Kinematics` / `Propagator` data
  types (`family.rs`), the registered symbol set including `dot = symbol!(...; Symmetric, Linear)`
  (`symbols.rs`), and `amplitude(&IntegralFamily) -> Atom` folding `Σ coeff · master_symbol`
  (`amplitude.rs`). Built on Symbolica 2.1.0.

---

## Notes on scope and speed

- **Scope:** one-loop closed-form reduction only — replicating the Denner / FJT / Tarasov recursions in
  a CAS and verifying numerically. Not a Laporta engine; general Laporta and ≥2-loop work is a separate
  future effort.
- **Speed (context, not a claim of advantage):** symbolic reduction is sub-millisecond across topologies
  — triangle scalar 0.016 ms / rank-1 0.097 ms / rank-2 0.272 ms; box 0.029 / 0.142 / 0.543;
  pentagon scalar 0.504 / rank-1 0.727 (generic massive kinematics). This is exact-rational symbolic
  algebra producing a formula in `d`, done once per topology + numerator structure. MadLoop's warmed-up
  per-point numerical eval is 13.6 µs; oneloop's value is analytic closed-form output, not raw
  per-point numerical speed. See [benchmarks](06-benchmarks.md).

## Known gaps

- **Massless internal propagators** (`m² = 0`, IR-divergent masters): the `reg_delta` limit interplays
  with the master's own IR structure; guarded (clean error), scoped as follow-up.
- **Exactly-singular Gram box** at fully rank-deficient kinematics: clean error rather than a value.
- **Higher-rank numerators (rank ≥ 3)** beyond triangle/box, and the standalone EVALUATE phase
  (wiring OneLOopBridge to put numbers on the masters), are open.

See [the frontier](04-frontier.md) for the current edge and [the app](05-app.md) for integration status.

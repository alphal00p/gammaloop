# Validation of the bmhvkit toolkit

Saved by the coordinator from the result returned by the validation agent (key toolkit_validator).
Run the checks with: `cd SCRATCHPAD/toolkit/validation && python run_validation.py`.

## Summary

I found no bug in bmhvkit, so the toolkit code is unchanged. The author's suite passes (31/31, 21 s). My own suite also passes: 18 checks in about 85 s, run with `cd SCRATCHPAD/toolkit/validation && python run_validation.py`.

REPORT.md was not written. The harness blocks subagents from writing report files, so its content is in key_results, errors_found_in_setup and open_issues. It needs to be saved by the parent or by someone with write permission.

How the checks work:
- **Independent reference.** I wrote my own matrix evaluator (vcommon.word_matrix). It builds every input string directly from the explicit gamma matrices and g5, without going through bmhvkit's numeric_value.
- **Algebra and traces.** These are compared with that evaluator at d = 6, 8, 10.
- **Tensor reduction (up to rank 6).** This is compared with exact averages over the unit sphere, computed monomial by monomial.
- **uv_pole.** This is compared with a separate Feynman-parameter pole extraction for 2-, 3- and 4-point integrands with masses, open strings and traces. I also checked:
  - squared propagators against mass derivatives;
  - that shifting the loop momentum leaves the result unchanged;
  - the bubble_exact rank-1 and rank-2 coefficients against values I derived by hand;
  - the omega option.

Literature check, chiral Yang-Mills in Feynman gauge (LP = arXiv:2610.38618, BMIS = arXiv:2004.14398). I worked out the convention map myself: our Delta_1 at c = R is identical to LP's breaking term, and LP's momentum-space functions equal our diagram sum divided by i. Every result agrees exactly, and both chiralities give mirror images:
- (i) the fermion self-energy pole;
- (ii) the fermion loop in the gauge two-point function, both the four-dimensional pole and the evanescent -1/3 S_R phat^2 gbar pole;
- (ii) the c G breaking, pole and finite part, including LP's evanescent finite piece;
- (iii) c psibar psi, from three diagrams. The gauge exchange between the fermion legs has no pole and gives the C_R - C_G/2 finite part. The two ghost-gauge diagrams give the whole pole (C_G/2 times the tree breaking vertex) and +C_G/4. The total, C_R - C_G/4, is LP's (5+xi)/6 C_R - xi C_G/4 at xi = 1;
- (iv) c G G: the f^{abc} part matches LP, and the anomalous d^{abc} eps part matches BMIS (6.18).

I confirmed the toolkit author's corrections to 03_setup.typ and found two new problems there:
- a Grassmann sign error in the Higgs-conjugate term of eq:ST;
- a conflicting definition of omega(G) in sec:integrals.

## Key results

- No bug found in bmhvkit, so the code is unchanged. The author's suite gives 31/31 pass; the validation suite (validation/run_validation.py) passes all 18 checks in about 85 s.
- Algebra: 200 random strings of 6-12 gammas (f/b/h kinds, up to 3 contracted pairs, free and four-dimensional indices, d-dim and 4-dim momenta, g5/PL/PR anywhere) agree with explicit matrices at d = 6, 8, 10, both simplify() and tr(). Also passing: 40 g5 traces of up to 8 gammas with contracted hatted pairs; tr(gbar^6 g5) and tr(gbar^8 g5) on 300 component sets each; 50 strings multiplied by eps tensors and metrics (eps-gamma, eps-eps and metric-eps contractions) at d = 6, 8; Chisholm reduction.
- tensor_reduce up to rank 6 (bar/hat/full scalar products, restricted k^2, slashed k) equals the exact average over the unit sphere, computed monomial by monomial, at d = 6, 8, 10.
- uv_pole equals an independent Feynman-parameter pole extraction (exact Gamma formula, exact simplex integrals) for 40 random 2-, 3- and 4-point integrands with masses, open strings and traces; the full eta-tagged pole N0 + eta N1 + eta^2 Nrest is compared. Also passing: shifting the loop momentum leaves N0 and N1 unchanged; squared propagators equal mass derivatives; bubble_exact rank-1 and rank-2 coefficients (pole, finite, log) match hand-derived values; the bubble_exact pole equals uv_pole.
- Convention map: our Delta_1 and eq:V-delta at c = R are identical to LP's breaking term and its rule (same field order in s psibar, g c <-> omega). LP's momentum-space function equals our diagram sum divided by i, which is (N0/eps + N1)/(16 pi^2). BMIS's p1 in Gamma_Gc is the gauge-field momentum.
- (i) Fermion self-energy: Gamma_div = C_R pbar-slash P_c / eps, with no evanescent pole. Agrees with LP and BMIS (5.10).
- (ii) Fermion loop in the gauge two-point function: Gamma_div = S_R [(2/3)(pbar pbar - pbar^2 gbar) - (1/3) phat^2 gbar] / eps. Agrees with LP and BMIS (5.12).
- (ii) Delta-inserted c A, with p the ghost momentum: (i S_R/3) pbar^nu (phat^2/eps + pbar^2 + phat^2). Agrees with LP eqn:Gc_viol_1_loop, including the evanescent finite pbar phat^2 term; the physical part of the exact bubble equals P N1.
- (iii) Delta-inserted c psibar psi, Feynman gauge, three diagrams. (a) Gauge exchange between the psi and psibar legs: group factor C_R - C_G/2, no pole, P N1 = i(pbar - pbar')slash P_c. (b), (c) External ghost -> ghost-gauge vertex -> internal ghost -> Delta vertex, with the gauge boson ending on the psi or the psibar leg: group factors +-(i/2) C_G; together they give the whole pole and +(i/4) C_G (pbar - pbar')slash P_c. Total: pole = (C_G/2)/eps times the tree Delta vertex (purely evanescent); finite = i(C_R - C_G/4)(pbar - pbar')slash P_c. This equals LP's (5+xi)/6 C_R - xi C_G/4 at xi = 1 and BMIS (6.10d), (6.24). The open ghost line gives no Grassmann sign.
- (iv) Delta-inserted c A A. The f^{abc} part, -(S_R/3) f^{abc} [(phat1^2 - phat2^2) gbar/eps - (pbar1^2 - pbar2^2) gbar + 2 pbar1 pbar1 - 2 pbar2 pbar2], agrees with LP GGc. The d^{abc} part has no pole and a finite (2/3) d_R eps^{mu nu rho sigma} pbar1_rho pbar2_sigma for c = R, which agrees with BMIS (6.18); the sign flips for c = L.
- Usage notes, not bugs: numeric_value returns a number when no Dirac string survives; Expr / int is not defined (use * sympy.Rational); with the explicit omega= option, omega must be the mass dimension of the integral. No HepKIT/Symbolica/gammaloop tool was used, so there is no gaps file for this key.

## Errors found in the setup section (all corrected in 03_setup.typ on 2026-11-27)

- eq:ST, the term (dF/dH^dagger)(dF/dsigma). The setup says derivatives act from the left by default, but F contains (s H^dagger) sigma with both factors Grassmann odd, so the left derivative with respect to sigma is -(s H^dagger). This term then gives -(s H^dagger) dS0/dH^dagger, and eq:ST-S0 fails for the H^dagger part. Corrected: (dF/dH^dagger)(delta_R F/delta sigma), or equivalently -(dF/dH^dagger)(delta_L F/delta sigma). Verified by hand; the other source terms are fine (rho and sigma^dagger stand to the left, zeta and R are bosonic).
- sec:integrals, the definition of omega(G) in the expansion operator. 'Superficial degree of divergence' read as k-power counting (deg_k + 4 - 2N) conflicts with counting the p_i and m_j in the numerator as degree one. Example: int pslash kslash/(k^2 (k+p)^2) has k-power omega = 1; with the literal rule only an odd term survives and the result is 0, but the true pole is -p^2/(2 eps). Corrected: omega(G) is the mass dimension of the integral, counting numerator momenta and masses (equivalently, keep every term whose vacuum integral has non-negative degree). Verified with validation/test_uvpole.py (test_omega_option); uv_pole does this by default.
- Text after eq:logpole, pole for b = a (confirming the toolkit author): the coefficient is a(a+1) M^4/(2 eps), not (a+1)(a+2) M^4/(2 eps). Check: int k^2/(k^2-M^2) = M^2 int 1/(k^2-M^2) gives M^4/eps for a = 1.
- Sentence after eq:matrix-fields (confirming the toolkit author): s wm = +i wm wm, not -i wm wm. Derived from s omega^a = -1/2 fg^{abc} omega^b omega^c and [Tg^b, Tg^c] = i fg^{bca} Tg^a.
- eq:gf-exact (confirming the toolkit author): the source part must read + sum_f (Rbar_f Psi_f - Psibar_f R_f) - sigma^dagger H + H^dagger sigma to reproduce eq:Lext with a left graded derivation (R is bosonic, sigma is odd). Checked term by term.

## Open issues / not validated

- REPORT.md in validation/ was not written: the harness blocks subagents from writing report files. Its content (test list, convention map, diagrams of (iii), setup errors, what remains unvalidated) is in this output and needs to be saved by the parent.
- Not validated against the literature: the Delta-inserted 4-point functions (psibar psi A c, Rbar psi omega omega, psibar R omega omega, c A A A, c A A A A) and all scalar/Yukawa sectors. The engine parts they rely on are validated generically.
- The group module was checked only by hand on tab:invariants and eq:SLR, which agree. The Yukawa, projection and rigid-invariance helpers rest on the author's tests.
- bubble_exact finite parts are checked only for rank <= 2 coefficients, through the author's Ward identity, and through P(finite) = P N1 for the Delta bubble.
- No general-xi test: the engine works with Feynman-gauge propagators only.
- N1 depends on the bookkeeping rule for vector-like and ghost-number-0 diagrams (agreeing with the author). It is unambiguous only after projection for Delta-inserted functions. Local finite pieces of ghost-number-0 functions, such as -(1/3) S_R pbar^2 gbar, are outside N1.


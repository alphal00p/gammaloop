# bmhvkit — checking toolkit for the one-loop BMHV renormalisation of the unbroken SM

Conventions: `report/sections/03_setup.typ` and `report/macros.typ`. Where the draft has errors, the toolkit uses the
corrected statements listed under "Deviations from the setup draft" below.

```
python=<scratch>/venv/bin/python
cd SCRATCHPAD/toolkit
$python tests/run_tests.py     # 31 tests, about 20 s  (test_* functions are pytest-compatible; pytest is not in the venv)
$python examples.py            # the examples of this README, about 2 s
```
To use it from a script: `sys.path.insert(0, ".../SCRATCHPAD/toolkit")`, then `import bmhvkit`. It needs only sympy
and numpy; Symbolica is not imported (see "Why not HepKIT/Idenso" below).

## Modules

| module | content |
|---|---|
| `bmhvkit.numeric` | back-end (i): explicit matrices in even integer d (4, 6, 8, 10, 12), exact interpolation in d |
| `bmhvkit.symbolic` | back-end (ii): symbolic strings of gbar, ghat, g5, PL, PR; contractions, traces, eps, the `eta` tag |
| `bmhvkit.loop` | one-loop UV engine: `T^(omega)` with auxiliary mass M, tensor reduction, vacuum integrals; returns N0, N1; exact massless bubble |
| `bmhvkit.group` | SM generators, structure constants, invariants, anomaly traces, Yukawa couplings, projection on invariant tensors |

## 1. Symbolic algebra (`bmhvkit.symbolic`, re-exported by `bmhvkit`)

Objects:
- `Idx("mu")` is a Lorentz index; `Idx("mu", four=True)` takes only 4-dimensional values.
- `Mom("p")` is a d-dimensional momentum; `Mom("p", four=True)` has a zero hat part.
- Linear combinations of momenta (`k + p`, `p1 - 2*p2`) are accepted wherever a momentum is.
- Slot kinds: `"f"` (d-dim, the default), `"b"` (bar), `"h"` (hat).

Constructors (each returns an `Expr`):
- gamma matrices: `G(x, kind)` (x an index or a momentum: `G(p)` is p-slash), `Gb(x)`, `Gh(x)`, `slash(v, kind)`;
- chirality: `g5()`, `PL()`, `PR()`, `proj("L"|"R")`, `one()`;
- tensors: `metric(x, y, kind)` gives g^{xy}, a vector component or a scalar product, according to whether x and y are indices or momenta; `vec(p, mu, kind)`, `dot(p, q, kind)`, `eps(a, b, c, e)` (4-dim, slots may be momenta).

`Expr` supports `+ - *` (the matrix product, in the order written) and multiplication by sympy scalars. Methods:
- `.simplify()` returns the canonical form;
- `.tr()` returns the Dirac trace (tr 1 = 4);
- `.phys()` and `.evan()` are the projections P and E of eq:Pproj (taken after simplification);
- `.subs(...)`, `.coeff(structure)`, `.collect()`, `.is_zero()`.

Utilities:
- `chiral_parts(e)` rewrites sum X(a + b g5) as {X: (c_L, c_R)}.
- `chisholm(e)` reduces runs of three or more barred gammas in open chains.
- `numeric_equal(a, b, dims, free_indices, trace_a, trace_b)` compares two expressions with explicit matrices, using random integer momenta and all components of the free indices.
- `numeric_value(e, Algebra(d), vectors, indices)` evaluates an expression with explicit matrices.
- `fresh()` returns a new dummy index.

Rules:
- Repeated indices are summed. Avoid clashes yourself when multiplying expressions; `fresh()` gives unique names. Names starting with `_` are reserved for the engine.
- Traces of metrics: g^mu_mu = `d`, gbar^mu_mu = 4, ghat^mu_mu = `-2*eta` (eq:eta). `d` and `eta` are sympy symbols exported by the package.

Canonical form produced by `.simplify()`:
- All contractions are done. The only exception is an index shared by an epsilon and a gamma matrix, which is kept.
- Every d-dim object is split into bar + hat, so only kinds b and h remain.
- There is at most one epsilon per term, because eps·eps is reduced with eq:epseps.
- In the Dirac string, all gbar stand left of all ghat, each group is sorted, and a single g5 sits at the right end. Projectors are expanded as (1 ± g5)/2; use `chiral_parts` to read the result in PL/PR form.
- The form is unique up to 4-dim Schouten identities, and up to the Chisholm identity for open chains with three or more barred gammas. For such cases use `numeric_equal`.

**Order of operations, which fixes the eta bookkeeping:**
1. A contracted pair that is not separated by g5 is contracted first, in d dimensions:
   `G(a)*X*G(a)` gives `d`, e.g. g^a gbar^4 g_a → (d-8) gbar^4.
2. A pair separated by g5 is split before g5 is moved:
   `G(a)*g5()*G(a)` gives `(-4 - 2*eta)*g5` (= (d-8) g5 at eta = eps).
3. d-dim momenta and metrics are split into bar + hat after the contractions.

```python
from bmhvkit import *
mu, nu, a = Idx("mu"), Idx("nu"), Idx("a"); k, p = Mom("k"), Mom("p")
(Gh(mu)*Gb(nu)*Gh(mu)).simplify()                   # (2*eta)*[gb(nu)]            i.e. -(d-4) gbar^nu
(G(a)*Gb(mu)*Gb(nu)*Gb(Idx("r"))*Gb(Idx("s"))*G(a)*g5()).tr()   # (-4*I*(d - 8))*eps(mu,nu,r,s)
(PR()*slash(k)*PL()).simplify()                     # (1/2)[gb(k)] - (1/2)[gb(k) g5]  = kbar-slash PL
e = G(mu)*slash(p)*G(nu)*slash(k)*PR()
numeric_equal(e.tr(), e, dims=(4, 6), free_indices=["mu", "nu"], trace_b=True)   # True
```
How results print:
- `pb(mu)` is pbar^mu, `ph.ph` is phat^2, `gb(mu,nu)` is gbar^{mu nu};
- a bracket `[gb(p) gh(q) g5]` is the Dirac string pbar-slash qhat-slash g5.

## 2. Explicit matrices (`bmhvkit.numeric`)

`Algebra(d)` provides:
- the matrices and metric: `.g[m]` (gamma^m), `.g5` (= i g^0 g^1 g^2 g^3, built from the first four matrices), `.PL`, `.PR`, `.one`, `.eta`;
- trace and epsilon: `.tr(M)` (normalised to tr 1 = 4), `.eps_up(a,b,c,e)` (eps^{0123} = +1);
- projections on a kind (`"b"`, `"h"` or `"f"`): `.gamma(m, kind)`, `.slash(v, kind)`, `.metric(kind)`, `.project(v, kind)`, `.project_matrix(M, kind)`, `.dot(u, v, kind)`.

Functions:
- `evaluate(alg, factors, fixed, trace)` evaluates a product of `("g", label, kind)`, `"g5"`, `"PL"`, `"PR"` and matrices, summing over repeated labels.
- `interpolate_in_d(func, dims)` returns the exact polynomial in d. Matrix entries are 0, ±1, ±i and momenta are integers, so the values are exact dyadic rationals.
- `exact(z)` converts a value to an exact sympy number.

## 3. One-loop UV engine (`bmhvkit.loop`)

```python
from bmhvkit.loop import Den, uv_pole, bubble_exact, tensor_reduce, substitute_momentum, ep, M
r = uv_pole(numerator, k, [Den(shift, mass2, power), ...], omega=None)   # -> UVResult
r.check()          # asserts that the M-dependence cancels
r.N0, r.N1         # T^(omega) G = i/(16 pi^2) (N0/eps + eta N1/eps + eta^2 Nrest/eps) + finite   (eq:N0N1)
```

Integrand:
- G = ∫_k numerator / ∏ ((k+shift)^2 - mass2)^power, with ∫_k = mu^{2eps} ∫ d^dk/(2pi)^d.
- The numerator contains every i, sign and coupling from the Feynman rules, and the -1 of a closed loop.
- The numerator may be an open Dirac string or a traced scalar. External momenta are `Mom` (d-dim, or `four=True`); external indices are `Idx` (d-dim, or `four=True`).

Steps:
1. Expansion: 1/((k+q)^2 - m^2)^n = Σ_j binom(-n,j) B^j/(k^2-M^2)^{n+j}, with B = 2k·q + q^2 - m^2 + M^2.
2. Terms are kept if their UV degree deg_k - 2·(power of 1/(k^2-M^2)) + 4 is ≥ 0. For homogeneous integrands this is T^(omega(G)). With `omega=n` and `mass_degrees={m: 1}`, the small-quantity degree ≤ n is used instead; pass masses as `m**2`.
3. Canonical split, then tensor reduction with the d-dim metric (eq:tensor-red). Only the denominators' k^2 stays unsplit: a numerator k^2 becomes kbar^2 + khat^2 → (4 - 2 eta) k^2/d.
4. Vacuum poles: the pole of ∫ (k^2)^a/(k^2-M^2)^b is i/(16pi^2) c_ab M^{2m}/eps, with m = a+2-b and c_ab = (-1)^{a+b+m}(a+1)!/(m!(b-1)!).
5. The numerator coefficient is taken at d = 4 with eta kept. The result is split by powers of eta into N0, N1 and Nrest; M-dependent leftovers go to `Mrem`.

Exact bubble:
- `bubble_exact(num, k, p, eta_mode="eps")` computes ∫_k num/(k^2 (k+p)^2) exactly in d (Feynman parameters) to O(eps^0).
- It returns `pole`, `finite` and `log`, with G = i/(16pi^2)[pole/eps + finite + log·L], L = log(-p^2/mu^2) (MS-bar).
- It is massless and two-point only. Use it to check finite parts and Ward/ST identities exactly.

Other functions: `tensor_reduce(e, k)` and `substitute_momentum(e, k, v)` (loop-momentum routing, k → k+p) are also exported.

```python
I = sympy.I; S = lambda x: I*slash(x)
# chiral fermion self-energy, Feynman gauge, unit generator (c = L)
num = (I*Gb(mu)*PL()) * (I*slash(p+k)) * (I*Gb(mu)*PL()) * (-I)
uv_pole(num, k, [Den(0,0), Den(p,0)]).check().N0        # chiral_parts -> {'[gb(p)]': (1, 0)} = pbar-slash PL ; N1 = 0
# fermion loop in <A A> and the diagram with the breaking vertex eq:V-delta (psi in k, out k+p)
NAA = -((I*Gb(mu)*PL()) * S(k) * (I*Gb(nu)*PL()) * S(k+p)).tr()
ND  = -((slash(k,"h")*PR() - slash(k+p,"h")*PL()) * S(k) * (I*Gb(nu)*PL()) * S(k+p)).tr()
uv_pole(NAA, k, [Den(0,0), Den(p,0)]).N0   # 2/3 (pb pb - pb.pb gb) - 1/3 ph.ph gb(mu,nu)
uv_pole(ND,  k, [Den(0,0), Den(p,0)]).N0   # (I/3) pb(nu) ph.ph
uv_pole(ND,  k, [Den(0,0), Den(p,0)]).N1   # (I/3) pb.pb pb(nu) + (I/3) pb(nu) ph.ph
bubble_exact(NAA, k, p).finite.phys()      # -13/9 pb.pb gb + 10/9 pb pb   (non-transverse: -1/3 pb.pb gb)
```
These results are self-consistent:
- p_mu D_AA^{mu nu} = i D_Delta^nu holds exactly (pole, finite and log);
- pbar_mu P[finite_AA] = i P[N1_Delta] (test `test_C_ward_identity_with_delta_vertex`).

Run time: a triangle with one Delta vertex takes about 0.1 s; the fermion box (8 gamma matrices plus projectors, trace of 703 terms) takes about 1.5 s.

## 4. SM group helper (`bmhvkit.group`)

Adjoint index a = 0..11 (G^1..G^8, W^1..W^3, B). Multiplet components:
- Q: 2·alpha + i;
- U, D: alpha;
- L, H: i;
- E: a single component.

Couplings are the symbols `gs, g, gp`; Yukawas are `yu, yd, ye`. Pass `couplings="unit"` to set the gauge couplings to 1.

- **Generators and structure constants.** `generators(f)` returns the twelve Tg^a_f (eq:Tg-explicit). `fg()` returns fg^{abc} (12×12×12). `su3_f()` and `su3_d()` give the SU(3) f and d.
- **Group invariants.** `Tr(*mats)`, `S_index(f, "G"|"W"|"B")`, `casimir(f)`, `casimir_adjoint(factor)`, `anomaly_D(a,b,c)` and `check_lie_algebra(f)`.
- **Yukawa couplings.** `yukawa_terms()` returns L_Yuk as entries {fbar, f, scalar ('H' or 'Hdag'), proj, coef[r][s][j]} (eq:LYuk); the vertex is i·coef·P. `yukawa_matrix(fbar, f, scalar, j)` returns one coefficient matrix. `check_yukawa_invariance()` checks rigid invariance.
- **Projection on invariant tensors.**
  - `project(X, basis)` gives the exact coefficients and the residual.
  - `adjoint_rank2_basis()` gives δ per factor.
  - `adjoint_rank3_basis()` gives fG, dG, epsW, the δδ_B placements and BBB.
  - `is_invariant(X, legs)` checks rigid invariance. Use `("adj",)` for an adjoint leg, `(f, +1)` for an index that transforms like the components of f, and `(f, -1)` for the conjugate. For the coefficient of a field monomial, use the conjugate of each field: fbar_r → `(f,+1)`, f_s → `(f,-1)`.

## Tests (tests/run_tests.py)

**A. Explicit matrices**, d = 4..10:
- Clifford relations, g5^2 = 1, g5 hermitian, {gbar, g5} = 0, [ghat, g5] = 0, {gamma, g5} = 2 ghat g5;
- hermiticity;
- eq:g5def;
- tr 1 = 4, tr g5 = 0, eq:tr5 and odd traces;
- eq:proj and eq:sandwich;
- eq:metric, eq:contr1 and eq:contr2 as exact polynomials in d;
- eq:epseps;
- factorisation tr(B H) = tr(B) tr(H)/4.

**B. Symbolic engine**:
- the same identities, checked symbolically;
- eq:phys-prop and the P/E projections;
- 25 random strings (f/b/h, indices, d-dim and 4-dim momenta, g5/PL/PR), canonical forms and traces compared with matrices at d = 4, 6, 8;
- g5-odd traces of 6 barred gammas with all components;
- a trace with full contractions compared with interpolation in d;
- Chisholm reduction.

**C. Loop engine**:
- vacuum-integral pole coefficients against the exact Gamma-function formula;
- tensor reduction, including the bar/hat projections;
- massive scalar integrals and the cancellation of M;
- QED: self-energy N0 = p-slash, vacuum polarisation (4/3)(pp - p²g), vertex gamma^mu;
- QCD: gluon plus ghost loop (5/3)(p² g - pp) in Feynman gauge, and no pole from the seagull;
- chiral self-consistency tests derived by hand in the docstrings:
  - chiral self-energy pbar-slash P_c with N1 = 0;
  - fermion loop in <AA>: evanescent non-transverse pole -1/3 phat² gbar, and non-transverse finite part -1/3 pbar² gbar from the exact bubble;
  - Ward/ST identity with the Delta vertex: exact in d, at the level of poles, and N1;
- exact B0 and scaleless integrals;
- the `omega` option.

**D. Group helper**:
- Lie algebra on every multiplet;
- every entry of tab:invariants, C_ad and eq:SLR;
- D^{abc} = 0 for all a, b, c;
- Yukawa invariance, together with a deliberately wrong assignment that must fail;
- projection of fg onto gs·fG + g·epsW, and invariance checks.

## Deviations from the setup draft (03_setup.typ)

The toolkit implements the corrected statements:
1. **Pole of the vacuum integral with b = a** (text after eq:logpole). The correct coefficient is i/(16pi^2) · a(a+1) M^4/(2 eps), not (a+1)(a+2)M^4/(2 eps). Check: ∫k²/(k²-M²) = M² ∫1/(k²-M²) gives M⁴/eps for a = 1. This is test `test_C_vacuum_integral_poles`.
2. **The BRST variation of the matrix ghost** (sentence after eq:matrix-fields). From sop gh^a = -½ fg^{abc} gh^b gh^c it follows that sop wm = +i wm wm, not -i wm wm. The same result follows from sop² Psi = 0 with sop Psi = i wm Pc Psi.
3. **eq:gf-exact.** With the stated left graded derivation, the source terms must be + Σ_f (Rbar_f Psi_f - Psibar_f R_f) - σ^† H + H^† σ to reproduce eq:Lext. The draft has + Psibar R and - H^† σ.
4. **The eta bookkeeping (eq:eta, eq:tensor-red-split) is not unique.** kbar² + khat² → (4 - 2 eta) k²/d, while an unsplit k² stays k². Similarly, g^a g5 g_a is -4 - 2eta or -d - 4eta depending on how g5 is moved. N1 is therefore well defined only once a splitting rule is fixed. The engine uses one fixed rule (see §1 and §3). With that rule, N1 ≠ 0 for purely d-dimensional (vector-like) diagrams, e.g. QED vacuum polarisation N1 = 2/3 p² g; this is irrelevant, because only poles are used there. For the Delta-insertion functions every hat comes from the BMHV projections, and P N1 is unambiguous. It agrees with the exact finite part (test C).
5. **Finite parts of ghost-number-0 functions.** The non-transverse finite part of the fermion loop in <AA>, -1/3 pbar² gbar per unit S_f, is not an eta/eps term in the setup bookkeeping (N1 = 0). It comes from 1/(d-2) in the full integral, and it is recovered through the Delta insertion: i P N1_Delta = pbar_mu P finite_AA.

## Not covered

- Diagram generation, symmetry factors and the Grassmann signs of ghost and source legs. The user writes integrands from the Feynman rules. HepKIT can enumerate diagrams; see `notes/code_hepkit_model_and_generation.md`.
- Group factors are not attached automatically to the Lorentz/Dirac part. Compute them with `group` and multiply them in.
- The exact (finite-part) integrator covers only massless bubbles. UV poles are fully general, including masses via `Den(…, m**2)` and higher powers.
- Open chains are canonical only up to 4-dim Schouten/Chisholm identities. Compare with `numeric_equal`, or use `chisholm`.
- Higher loops: there is no subdivergence subtraction (R-bar) and no IR rearrangement beyond the single-M T-operator.
- Performance is that of plain Python and sympy. Traces with more than about 10 barred gammas plus g5 are slow (the epsilon recursion produces non-minimal forms).
- Rigid-invariant bases are provided only for rank-2 and rank-3 adjoint tensors. Other bases must be passed to `project`.

## Why not HepKIT/Idenso for the symbolic engine

The Symbolica tensor module leaves mixed mink(4)/mink(D-4) words untouched. For example, `tr(gbar^mu ghat^nu)` and `tr(g5 ghat ghat)` come back unchanged. This is recorded in `notes/gaps/toolkit_author.md`. The BMHV algebra is therefore implemented directly; it is cross-checked against explicit matrices.

# The on-shell massless-leg frontier

`oneloop` is a *massive-configuration* reducer: its per-topology recursions are
built on mass-derivative / Tarasov-style relations that divide by Gram and
(modified) Cayley determinants. That machinery is exact and complete on generic
kinematics (see [the reduction algorithm](02-reduction.md)), but it has one
structural weakness — it degrades when external invariants (or internal masses)
vanish. That regime, **on-shell massless external legs**, is not exotic: it is
*universal* in real collider processes, where external gluons, photons, and light
quarks all satisfy `p² = 0`. This document maps that frontier precisely, describes
the fix (`reg_delta` off-shell regularization), gives the validated numbers, and
states honestly what remains a genuine gap.

## Why the massive-configuration method breaks on-shell

The reducer expresses tensor-integral coefficients through inverse Gram / Cayley
matrices of the external kinematics. For a triangle, the modified Cayley matrix is

```
Y = [[2m0², s01−m0²−m1², s02−m0²−m2²],
     [ ...            , 2m1², s12−m1²−m2²],
     [ ...            , ...             , 2m2²]]
```

(schematically; `s_ij` are the pairwise external invariants, `m_i²` the internal
masses). When two external legs go on-shell and massless, a full row/column of the
associated Gram/Cayley structure collapses and `det(Y) → 0`. The recursions then
try to divide by that vanishing determinant. This is the classic
**inverse-Gram-determinant / exceptional-kinematics** problem.

Concretely, on the massless on-shell triangle `C0(0,0,s;0,0,0)` the modified
Cayley (`Y_ij = m_i² + m_j² − s_ij`, here all masses 0) is
`Y = [[0,−s01,−s02],[−s01,0,−s12],[−s02,−s12,0]]` with lex invariants
`[s01,s02,s12] = [0,0,100]`; it has a zero row (the two on-shell legs) and
`det(Y)=0` — degenerate.
The same mechanism strikes *pinched sub-topologies*: on-shell massless external
legs make a pinched sub-bubble's external leg the null on-shell leg, so its `1×1`
sub-Gram is `q² = 0` (singular Gram), or a pinched bubble `B0(0;m,m)` has a
`2×2` Cayley with zero determinant (degenerate Cayley).

Before the fix, the reducer's guards fired *cleanly* here — panicking with
`"singular Gram matrix"` or `"degenerate Cayley reduction: exceptional
kinematics"` (`reduce.rs`) rather than returning garbage. The integrals are
perfectly well-defined; the old algorithm simply did not cover that corner and
refused rather than producing a wrong number.

## The frontier map

The behaviour is more nuanced — and more encouraging — than "on-shell breaks it."
`frontier_map.rs` scans a massive-internal (`m²=1`) triangle and box, varying the
number of on-shell (zero-invariant) legs and the numerator rank, with each case
wrapped in `catch_unwind` so one run reports the full OK/PANIC map. The picture
*before* the `reg_delta` fix:

| Topology | Kinematics | Rank | Result |
|---|---|---|---|
| Triangle | 0 on-shell legs (all off) | rank-2 `(k·q1)²` | OK |
| Triangle | 1 on-shell leg | rank-2 `(k·q1)²` | OK |
| Triangle | 2 on-shell legs | rank-2 `(k·q1)²` | PANIC (singular Gram) |
| Triangle | 3 on-shell legs | rank-2 `(k·q1)²` | PANIC |
| Triangle | 2 on-shell `[0,0,0.4]` | scalar | OK |
| Triangle | 2 on-shell `[0,0,0.4]` | rank-1 `dot(k,q1)` | OK |
| Triangle | 2 on-shell `[0,0,0.4]` | rank-1 `k²` | OK |
| Triangle | 2 on-shell `[0,0,0.4]` | rank-2 `(k·q1)²` | PANIC (singular Gram) |
| Triangle | 2 on-shell `[0,0,0.4]` | rank-2 `(k·q1)(k·q2)` | PANIC |
| Triangle | 2 on-shell `[0,0,0.4]` | raised power `[2,1,1]` | PANIC (degenerate Cayley) |
| Box | all off-shell | rank-2 `(k·q1)²` | OK |
| Box | 2 on-shell legs | rank-2 `(k·q1)²` | **OK** |

Reading the map:

- **A triangle tolerates exactly ONE on-shell leg.** It takes TWO on-shell legs to
  break rank-2. Off-shell, the triangle is non-degenerate (e.g. `det Y = -0.32`
  for the gg→h massive-top case) and everything reduces.
- **The boundary is rank-1 → rank-2.** At two on-shell legs the scalar and both
  rank-1 numerators (`dot(k,q1)`, `k²`) reduce; rank-2 and raised powers do not.
- **Boxes are robust.** The box reduces at two on-shell legs even at rank-2. Its
  richer topology pinches to *triangles*, not directly to null-leg bubbles, so it
  avoids the degenerate sub-Gram. The degenerate frontier is worst for the
  triangle — the topology with the fewest legs.

So the gap to real processes is narrow and specific: **the multi-on-shell-leg
triangle at rank ≥ 2 (or with a raised propagator power)** — exactly the
IR / exceptional-kinematics corner that needs special reduction.

## The fix: off-shell `reg_delta` regularization

The massive-configuration method works perfectly wherever the Gram/Cayley
determinants are non-zero, and — crucially — its coefficients are **exact
rationals**, so the `1/delta` inverse-Gram poles are explicit. That gives a clean
route around the degeneracy:

1. Replace each degeneracy-causing zero invariant / mass with a symbol
   `oneloop::reg_delta` (a small off-shellness `p² = delta`).
2. Call the ordinary `reduce_core()` — the massive reducer, unchanged.
3. Substitute `delta → 0` in the resulting coefficients *and* master arguments,
   dropping zero terms.

This is **exact**, not an approximation: the `1/delta` inverse-Gram poles cancel
algebraically in the sum `Σ cᵢ Mᵢ`. Symbolica's exact rational arithmetic keeps the
cancellation perfectly clean — there is no numerical instability to fight, so the
`delta → 0` limit is taken symbolically for free.

`reduce()` implements exactly this (see [the overview](01-overview.md) and
[the reduction algorithm](02-reduction.md)): if any `kinematics.invariants` entry
`is_zero()` it delegates to `reduce_regularized`, a thin wrapper that regularizes
the vanishing invariants with `reg_delta`, calls the unchanged core reducer, and
takes the `delta → 0` limit; otherwise generic kinematics take the fast path
unchanged. **Blast radius: zero core changes.** No modification
to `reduce_core`, `reduce_cayley`, `gram_solve`, or `degenerate_coeffs`. Degeneracy
is detected by a pre-check (`is_zero`) rather than `catch_unwind`, because a caught
`reduce()` panic poisons Symbolica's global state — so one must not catch a
`reduce` panic and then reduce again in the same process.

A symbolic probe confirms the mechanism: `reduce()` accepts a *symbolic* invariant
`delta` without panicking, and the coefficients come out polynomial in `delta`
(finite at `delta=0`). For the gg→h rank-2 `(k·q1)²` case:

```
1/4 · delta² · C0(delta, 2/5, delta) + (…delta…)·B0(delta) + (…)·B0(2/5)
```

which at `delta=0` collapses to `(1/20)[B0(0;1,1) − B0(2/5;1,1)]`. The exact
rational arithmetic cancels the `1/delta` poles automatically. The same
regularization also covers **massless internal lines** (regularize `m² = delta`);
multiple massless legs share a single `delta`.

## Validated cases

All numbers below come from the MadLoop reproduction record and the deep-dive
sweep; see [benchmarks](06-benchmarks.md) and `../benchmarks/README.md`.

### gg→h massive-top triangle, rank-2, on-shell — matches MadLoop to 3.4e-7

The headline validation. Massive top internal (`m²=1`), legs `(delta, delta, 0.4)`,
rank-2 numerator `(k·q1)²`.

**Off-shell ladder + extrapolation.** Reducing at small off-shellness `p²=delta`
and taking `delta → 0`, the finite (ε⁰) part:

| delta | finite part (ε⁰) |
|---|---|
| 1e-1 | +4.2e-4 |
| 1e-2 | −2.89e-3 |
| 1e-3 | −3.415e-3 |
| 1e-4 | −3.469e-3 |
| 1e-5 | −3.4742e-3 |

Fitting `A + B·delta·ln(delta) + C·delta` gives the on-shell limit
`A = −3.475174e-3`. There is no numerical breakdown down to `delta=1e-5` — the
exact-rational coefficients keep the cancellation clean. The UV pole
`ε⁻¹ = 0.25·delta → 0`, physically correct since the `(k·q1)²` UV divergence
tracks `q1² → 0`.

**Independent, convention-free cross-check.** On-shell the trace piece `~ q1² = 0`,
so the integral equals `C0 · ⟨(P·q1)²⟩_F` (F-weighted Feynman-parameter simplex
average), with `C0` from OneLOop: `−3.474834e-3`.

**Difference between the two:** the off-shell→0 extrapolation `−3.475174e-3` vs the
independent `−3.474834e-3` differ by **3.4e-7**. Match.

**Symbolic (limit-free) route.** Evaluating the symbolic-`delta` collapse
`(1/20)[B0(0;1,1) − B0(2/5;1,1)]` directly: `ε⁻¹ = 0` (UV pole cancels),
`ε⁰ = −3.4748337e-3`, vs the validated `−3.474834e-3` — difference **3.2e-10**.
The analytic route avoids even the numerical extrapolation residual.

### Deep-dive: 2 / 3 / all massless legs (`massless_legs.rs`)

The deep-dive sweep runs 13 on-shell configurations (plus 4 off-shell
path-independence checks), each in its own process (to avoid the `catch_unwind`
poison). **Zero panics** — the fix handles every case.
Results, all physically sensible and verified by identity/limit/hand-derivation:

- **2 massless legs, massive internal (gg→h)** — the baseline above. scalar
  `C0(0,s,0;1,1,1)`; rank-1 `½B0(0) − ½B0(2/5)`; rank-2
  `1/20·B0(0) − 1/20·B0(2/5) = −3.4748e-3`. Validated to 1e-10.

- **3 massless legs, massive internal** `C0(0,0,0;m,m,m)` (collinear point):
  scalar `C0(0,0,0;1,1,1)` (finite collinear master); rank-1 → 0; rank-2 → 0;
  `k²` → `C0(0,0,0;1,1,1) + B0(0;1,1)` (exact `k² = D0 + m²` identity). Rank ≥ 1
  external-momentum numerators vanish because in the fully-collinear limit all
  external momenta become proportional and shrink (`q_i·q_j → 0`), so the
  `⟨dot(k,q_i)⟩` projections vanish. **Path-independence confirmed:** reducing
  rank-2 at small *off-shell* invariants (direct path, no regularization) → 0 as
  `delta²`, the same for symmetric (`0.1 → −8.6e-4`, `0.01 → −8.4e-6`) and
  asymmetric (`0.01 → −1.2e-4`, `0.001 → −1.2e-6`) invariants — so the 0 is genuine,
  not a single-`delta` artifact.

- **2 massless legs, MASSLESS internal (e+e→dd̄, off-shell leg `s=1`):** scalar
  `C0(0,s,0;0,0,0)` (IR-divergent master, OneLOop handles); rank-1 `dot(k,q1)`
  → `½B0(0;0,0) − ½B0(s;0,0)` (hand-verified, `B0(0;0,0)=0`); **rank-1 `dot(k,q2)`
  → 0** (this previously *panicked*: both bubble pieces are scaleless, so the
  correct value is 0); rank-2 `1/8·B0(0;0,0) − 1/8·B0(s;0,0)`. The fix therefore
  also cures the massless-*internal* e+e→dd̄ case — `dot(k,q2)` and the tensor
  numerators now return correct values instead of panicking.

- **All massless (external + internal):** scalar → `C0(0,0,0;0,0,0)`, which is
  scaleless = 0 (OneLOop evaluates it to 0; the reducer keeps the master symbol
  rather than auto-simplifying); rank-2 → 0. In dim reg the answer is trivially 0 —
  no special reducer is needed here.

## The remaining gap

Honest scope. The regularization is powerful but not universal.

- **Continuity caveat.** The single-shared-`delta` diagonal limit gives the correct
  answer *when the integral is continuous at the on-shell point*. This holds for
  massive-internal finite integrals, and for the IR-regulated massless-internal case
  (where the ε-pole *coefficients* are continuous in the invariants). It can fail
  for **genuinely singular / threshold configurations** where different invariants
  must vanish at *different rates* — those need the Denner–Dittmaier systematic
  expansion. All 13 tested configs work; the rigorous general solution is
  DD / AMFlow-style.

- **Massless internal propagators / threshold.** Even with the delta fix cured for
  the tested e+e→dd̄ tensors, genuinely IR-divergent masters
  (`B0(0;0,0) = 1/ε_UV − 1/ε_IR`, `A0(0) = 0`, scaleless `C0`) still need the
  external master evaluator (OneLOop) to supply the correct dim-reg values;
  `oneloop` reduces to them but does not evaluate them (see
  [the app](05-app.md) and [the overview](01-overview.md)). Exactly-singular Gram
  boxes at all-massless kinematics remain guarded (clean-error, no garbage).

- **What `oneloop` does NOT need.** Fully-massless one-loop (external *and*
  internal) is scaleless = 0 in dim reg — no computation required. Tools like
  `vakint` / AMFlow target massive *vacuum* / multi-loop integrals (mass the only
  scale, no external momenta), a different regime; they are not needed for on-shell
  massless legs at one loop, which the `delta`-regularization already covers for
  2, 3, and the collinear/scaleless limits.

## Literature

The fix — *regularize the degeneracy with a small parameter, reduce, take the
limit* — is a well-known and rigorously developed idea:

- **Denner & Dittmaier, "Reduction schemes for one-loop tensor integrals,"
  Nucl. Phys. B734 (2006) 62, [hep-ph/0509141].** Exactly our case: they expand
  tensor coefficients about limits of vanishing Gram determinants (the
  inverse-Gram-determinant / exceptional-kinematics problem) to reduce to scalars.
  Their systematic Taylor expansion in the small Gram exists to avoid *numerical*
  instability; `oneloop` gets the same limit for free via Symbolica's exact rational
  arithmetic — same math, cleaner route. This is the rigorous general solution for
  the singular/threshold configs the single-`delta` limit does not cover.

- **Auxiliary Mass Flow (AMFlow), Liu–Ma–Wang, PRD 99 (2019) 071501;
  [arXiv:1711.09572] + package.** The rigorous, multi-loop realization of the
  "small mass → 0" intuition: introduce an auxiliary mass, compute in the simple
  (large-mass / vacuum) limit, then *flow* the mass back to physical via
  differential equations. The `reg_delta` external-leg regularization is the
  one-loop symbolic analog — regularizing external `p²` instead of internal masses.

- **Tarasov dimensional recurrence (`d → d+2`), hep-th/9606018.** The dimension-shift
  relations that underlie the mass-derivative recursions used here, and the natural
  route for a rigorous degenerate-kinematics reduction. Also the classic textbook
  precursor: mass regularization of IR divergences (give massless particles a small
  mass, take `m → 0`).

---

*See also:* [01-overview.md](01-overview.md) (scope and design),
[02-reduction.md](02-reduction.md) (the algorithm and `reduce()` wrapper),
[03-numerators.md](03-numerators.md) (tensor-numerator handling),
[06-benchmarks.md](06-benchmarks.md) and [../benchmarks/README.md](../benchmarks/README.md)
(full validation record), [CHANGELOG.md](CHANGELOG.md).

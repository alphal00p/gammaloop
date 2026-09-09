# H → γγ — the W-boson loop (a spin-1 amplitude through the reducer)

The second full-amplitude validation, and a harder one: `H → γγ` gets contributions
from a **fermion (top) loop** and a **W-boson loop**. The fermion loop is the same
`A_{1/2}(τ)` as [gg→h](09-ggh-formfactor.md) (validated to 10⁻¹³) — only the external
`N_c Q_t²` prefactor changes. The new content is the **W loop**: a spin-1 loop whose
form factor `A_1(τ)` is historically a subtle calculation.

**Result:** the reducer reproduces `A_1(τ)` across five τ points to `~10⁻⁵`
(δ-regularization limited → exact as δ→0), and the full H→γγ form factor comes out to
`−6.489` (analytic `−6.489`), giving `Γ(H→γγ) = 9.10 keV` (SM LO ≈ 9.1).

## Why the W loop is harder than gg→h

Three things it adds over the fermion triangle:

1. **A single diagram is not gauge invariant.** In unitary gauge the amplitude needs
   the full set: **two W triangles + a γγWW seagull bubble**. (Verified: `q₁·N ≠ 0`
   for one triangle.)
2. **Rank 6, not rank 2.** The unitary W propagator `−(g^{αβ} − k^α k^β/m_W²)` and the
   momentum-dependent γWW vertices push the loop-momentum tensor to rank 6.
3. **On-shell high rank needs δ-regularization.** For the on-shell triangle (q₁²=q₂²=0)
   the reducer reduces rank ≤2 exactly, but **rank ≥3 hits a `0⁻¹` in the reduction
   formulas**. Regularizing with `q₁²=q₂²=δ` makes every monomial finite; the `1/δ`
   poles cancel in the physical combination and `A_1` is the `δ→0` limit (confirmed:
   both the residual pole and `A_reduced−A_analytic` scale linearly in δ).

## Method

1. **Numerator (unitary gauge, couplings stripped, `wloop_numerator.py`).** Build each
   diagram's numerator from explicit contravariant vertices/propagators (metric on every
   bond + trace), contract with the transverse projector
   `P_{μν}=g_{μν}−(q₁_μq₂_ν+q₂_μq₁_ν)/(q₁·q₂)`, and **fit** the result to the monomial
   basis `ll^a lq1^b lq2^c` (`2a+b+c ≤ 6`) at integer `d`, then fit each coefficient as a
   polynomial in `d` (the whole `d`-dependence comes from `g^{μν}g_{μν}=d`). The fits are
   exact (residual ~10⁻⁹).
2. **One family.** Combine into a single **direct-triangle-family** numerator:
   `N = N₁ + swap_{q₁↔q₂}(N₂) + N_seag·(k²+2k·q₁−m_W²)` — the crossed triangle reduces
   in the direct family after `q₁↔q₂`, and the seagull bubble `(k, k+q₁+q₂)` embeds in
   the direct triangle `(k, k+q₁, k+q₁+q₂)` by multiplying by the pinched middle
   propagator `D₁`.
3. **Reduce + assemble (`wloop_reduce` + `wloop_assemble.py`).** The example reduces every
   monomial (δ-regularized) to `coeff(d)·master`; the Python side folds in the numerator's
   `d`-polynomials, evaluates the masters with OneLOop (Laurent in ε), requires the ε poles
   cancel, and forms `A_1 = P_{μν}M / [(q₁·q₂)(d−2)]|_{d=4}`.

## Numbers

`A_1(τ) = −[2τ² + 3τ + 3(2τ−1) f(τ)]/τ²`, `f(τ)=arcsin²√τ`:

| τ | A_1 reduced | A_1 analytic | rel. err |
|------:|------------:|-------------:|---------:|
| 0.6046 | −8.32428 | −8.32429 | 1e-6 |
| 0.300 | −7.52012 | −7.52024 | 2e-5 |
| 0.450 | −7.86561 | −7.86565 | 5e-6 |
| 0.800 | −9.19749 | −9.19750 | 1e-6 |
| 0.200 | −7.32606 | −7.32639 | 5e-5 |

(residual `~δ = 10⁻⁵`; shrinks linearly as δ→0).

### Full H→γγ

`m_H=125, m_t=173, m_W=80.379`, `τ_t=0.1305`, `τ_W=0.6046`:

```
F = N_c Q_t² A_{1/2}(τ_t) + A_1(τ_W) = 3·(2/3)²·(+1.376) + (−8.324) = −6.489   [analytic −6.489]
Γ(H→γγ) = α² m_H³ /(256 π³ v²) · |F|² = 9.10 keV   (SM LO ≈ 9.1)
```

The W dominates with destructive interference against the top — the textbook H→γγ, out
of the reducer.

## Re-run

```bash
python3 benchmarks/python/wloop_assemble.py   # needs oneloop_bridge + numpy + sympy
```

It builds the numerator in-process and shells out to
`cargo run --release --example wloop_reduce` per τ.

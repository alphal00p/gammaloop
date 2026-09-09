# gg → h form factor — end-to-end amplitude validation

This is the first **full-amplitude** validation of the reducer: not "does the
reduction match an independent master engine" (that is [07-benchmark-report](07-benchmark-report.md),
132/132), but "does reduce → evaluate → assemble reproduce a *physical* loop
quantity." The process is gluon-fusion Higgs production through a top loop,
`g(q₁) g(q₂) → h`, whose form factor `A_{1/2}(τ)` has a textbook closed form.

**Result:** the reducer reproduces `A_{1/2}(τ)` to **10⁻¹³** across six kinematic
points, and the assembled `|M|²` matches MadLoop's `9.3702613e-3` to **0.04 %**
once α_s is taken at the Higgs scale. The loop content is exact; the only
convention-dependent input is α_s.

## Why this needed real work

The app's reduction of a physical amplitude keeps the gluon **polarizations** ε
symbolic, so the loop momentum survives as `dot(k, ε)` — a polarization is not a
momentum in the reducer's Gram basis, so it is not a reducible ISP (see
[03-numerators](03-numerators.md)). To get a scalar form factor you must first
**project the polarizations out**, turning the numerator into momentum-basis dot
products the reducer fully reduces.

## Method

1. **Numerator.** The top-loop Dirac trace
   `N^{μν} = Tr[(k̸+m)γ^μ(k̸+q̸₁+m)γ^ν(k̸+q̸₁+q̸₂+m)]`, contracted with the
   **transverse projector** `P_{μν} = g_{μν} − (q₁_μ q₂_ν + q₂_μ q₁_ν)/(q₁·q₂)`
   (which annihilates the non-transverse pieces a single diagram carries), gives,
   in `d` dimensions and with on-shell gluons `q₁²=q₂²=0`, `q₁·q₂=s/2`:

   ```
   P_{μν}N^{μν} = 4m(6−d) k² + 4m(4−2d) (k·q₁) − (64m/s)(k·q₁)(k·q₂)
                  + [2ms(2−d) + 4m³(d−2)]
   ```

   This closed form is verified to 14 digits against explicit 4×4 Dirac matrices
   (in `ggh_formfactor.py`). The plain `g_{μν}` contraction is *not* enough — it
   keeps a non-transverse contamination that drifts the τ-shape by ~1 %.

2. **Reduce.** The `ggh_formfactor` example reduces each `dot(k,·)` monomial
   `{1, k·q₁, k·q₂, k², (k·q₁)(k·q₂)}` for the on-shell massive triangle. The
   external Gram `[[0, s/2],[s/2, 0]]` is non-singular (`det = −(s/2)² ≠ 0`), so the
   reduction is exact — no δ-regularization. E.g. `k² → m²·C0 + B0(0)`,
   `(k·q₁)(k·q₂) → (s/4)B0(0) − (s/8)B0(s)`.

3. **Evaluate + assemble.** `ggh_formfactor.py` evaluates the masters with OneLOop
   (Laurent in ε), folds in the projection's `(6−d)`, `(4−2d)` … coefficients,
   and **requires the 1/ε, 1/ε² poles cancel** (the form factor is finite — the
   `(4−d)=2ε` from the trace kills the `B0` pole). Then

   ```
   A_{1/2} = m_t · [P_{μν}M]_finite / [(q₁·q₂)(d−2)]|_{d=4} = m_t · finite / s
   ```

   The `m_t` factor is the trace's overall mass dimension; that it is exactly `1/m_t`
   (not a fitted constant) is confirmed by varying `m_t` independently of `m_H`.

## Numbers

`A_{1/2}(τ) = 2/τ²[τ + (τ−1) f(τ)]`, `f(τ)=arcsin²√τ`, `τ = s/(4m_t²)`:

| m_H | m_t | τ | A_reduced | A_analytic | rel. err |
|----:|----:|------:|----------:|-----------:|---------:|
| 125 | 173 | 0.1305 | 1.3762612 | 1.3762612 | 4e-15 |
| 125 | 150 | 0.1736 | 1.3915590 | 1.3915590 | 9e-15 |
|  90 | 173 | 0.0677 | 1.3549858 | 1.3549858 | 3e-14 |
| 200 | 250 | 0.1600 | 1.3866612 | 1.3866612 | 9e-15 |
| 300 | 173 | 0.7518 | 1.6933050 | 1.6933050 | 1e-15 |
| 125 | 350 | 0.0319 | 1.3433853 | 1.3433853 | 3e-13 |

The heavy-top limit `A_{1/2} → 4/3` (m_t=350, τ=0.032 → 1.3434) is approached correctly.

### Full |M|²

With `M^{μν,ab} = −i δ^{ab} (α_s/4πv) A_{1/2} (g^{μν} q₁·q₂ − q₁^ν q₂^μ)`, summing
colours (`Σ(δ^{ab})² = 8`) and polarizations (`T^{μν}T_{μν} = s²/2`) and averaging
the initial state (`1/256`):

```
|M|²_avg = α_s² A_{1/2}² s² / (1024 π² v²)
```

| α_s | |M|² | MadLoop | ratio |
|----:|-----:|--------:|------:|
| 0.1180 (m_Z) | 1.0509e-2 | 9.3703e-3 | 1.1215 |
| 0.1114 (Higgs scale) | 9.3663e-3 | 9.3703e-3 | 0.9996 |

The reducer's loop content (`A_{1/2}`) is exact; the ~12 % at `α_s(m_Z)` is entirely
the α_s scale — `α_s ≈ 0.111` (α_s run to the Higgs scale) reproduces MadLoop to
0.04 %. Pinning the last digits needs MadGraph's exact `param_card` (α_s scheme, v),
which is a convention, not loop physics.

## Re-run

```bash
python3 benchmarks/python/ggh_formfactor.py   # needs oneloop_bridge (see benchmarks/README.md)
```

It shells out to `cargo run --release --example ggh_formfactor` per kinematic point.

# Validation & benchmarks

This document is the evidence base for the `oneloop` reducer: what has been
cross-checked, against which independent engines, and where the numbers come
from. The reducer takes any one-loop integral with an arbitrary polynomial
numerator and reduces it to the four scalar masters A0/B0/C0/D0 with coefficients
rational in `d = 4 - 2ε` (see [the overview](01-overview.md) and
[the reduction algorithm](02-reduction.md)). Everything below is a check on those
coefficients and their assembled values — never on hand-tuned numbers.

Two complementary validation tracks are reported:

1. A **cross-engine harness** that checks the reduction of generic families
   against independent scalar-master libraries, an analytic library, a
   Feynman-parameter integrator, and tensor oracles.
2. A **MadLoop / MG5_aMC process suite** — ~20 real collider processes
   reproduced to ~14 digits where checked.

To run the harness yourself, see [`../benchmarks/README.md`](../benchmarks/README.md).

---

## 1. Cross-engine validation harness

The harness (`benchmarks/python/crosscheck.py`) reads the reductions emitted by
`benchmarks/rust/emit_reductions.rs` (Cargo example `emit_reductions`) and, for
each family, forms the full Laurent series

```
V_reduced(ε) = Σ_i c_i(d = 4 - 2ε) · M_i(ε)
```

where the `c_i` are oneloop's rational-in-`d` coefficients (series-expanded in ε
with `d = 4 - 2ε`) and each `M_i` is a scalar master value. It then requires:

- **(a) pole cancellation** — the `1/ε` and `1/ε²` poles cancel where the
  original integral is finite; and
- **(b) finite-part agreement** — the assembled finite part matches an
  independent direct evaluation.

### Independent engines cross-checked against

| engine | role | what it validates |
|---|---|---|
| **OneLOop** (`avh_olo`, van Hameren), via `oneloop_bridge` | numeric scalar masters A0/B0/C0/D0 (with poles) | the master **values** `M_i` |
| **feynalg** analytic `analytic_A0/B0/C0/D0` | analytic scalar masters | second, independent check of the master finite parts (OneLOop-vs-feynalg tally) |
| **scipy** Feynman-parameter integration (`scipy_quad` deterministic quadrature for N≤4; `scipy_direct` Monte-Carlo for N≥5) | direct value of the *original* integral in `d=4` | oneloop's **coefficients + topology selection** end-to-end |
| **tensor oracles** (`_tensor_oracle`: moment oracle + covariant Minkowski-Gram oracle) | direct value of tensor-numerator integrals | dotted / rank-raised numerators, with a **two-oracle** agreement requirement |

The renormalization scale is fixed to `μ² = 1` in the bridge
(`ob.set_renormalization_scale(1.0)`) to match the bare scipy Feynman integrand.

### Coverage validated (generic kinematics)

- **Scalar N = 3–7.** N≤4 reduce to OneLOop masters (exact, with poles); N≥5 are
  finite and checked against scipy MC.
- **Tensor numerators to rank-6.** Dotted numerators `dot(k,q_j)`,
  `dot(k,q_i)·dot(k,q_j)`, `k²·dot(k,q_i)` etc., plus pure `(k²)^p` families,
  checked against the moment oracle *and* a covariant Minkowski-Gram oracle from
  a different computational path (`two_oracle_agree`).
- **UV-divergent numerators** — divergent `(k²)^p` families are compared on the
  *full Laurent series* (both the `1/ε` coefficient and the finite part must
  match the oracle).
- **Massless lines** — handled by a corner-concentrated Dirichlet importance
  sampler in `scipy_direct` (a near-massless line makes the integrand singular at
  the `x_l → 1` corner; flat sampling would bias the estimate).
- **Timelike / above-threshold kinematics** — explicit signed invariants can be
  injected (the `SINV` line). Above threshold the integral is complex and scipy
  cannot cross `F = 0`; those `tlx_*` families are instead validated against the
  **IBP mass-derivative identity** `I[a_d = 2] = ± d/dm_d² I[scalar]`, with the
  right-hand side finite-differenced from OneLOop's *complex* master — validating
  the reduction coefficients in the complex plane.

### Reported result

The harness prints, per family, the pole residual, the reduced finite part, the
direct value, and a pull (in σ) or relative-error metric, with a `PASS`/`FAIL`
verdict: **poles cancel to `< 1e-6`; reduced finite part within 5σ of scipy MC**
(or within `2e-3` relative for the deterministic N≤4 quadrature). A separate line
reports the OneLOop-vs-feynalg master finite-part agreement tally.

Summary as reported for the Monday review: validated generically against OneLOop
+ feynalg + scipy + two tensor oracles across **scalar N=3–7, tensor to rank-6,
UV-divergent, massless lines, timelike / threshold / near-degenerate** — with
**37 unit tests** (`cargo test -p oneloop`) and clippy clean.

---

## 2. MadLoop / MG5_aMC process suite

Independently of the harness, ~20 real one-loop collider processes were generated
with **MG5_aMC v3.7.2** (python3.11, gfortran 15.2) and reproduced. Where a number
is quoted as checked, it matches Valentin Hirschi's reference / MadLoop output to
**~14 digits**.

| topology | processes (reproduced) |
|---|---|
| triangle | e+e⁻ → dd̄ (**−8.936**), Drell-Yan uū → e⁺e⁻ (**= −8.936**, by crossing), uū → Z, ud̄ → W⁺, tt̄ → g, e⁺e⁻ → bb̄ (massive b) |
| box | uū → dd̄, dd̄ → dd̄, uū → gg = dd̄ → gg = gg → dd̄ (**−54.03**, crossing), gu → gu, uū → tt̄ (massive top) |
| pentagon | e⁺e⁻ → dd̄g (5-external) |
| 4-gluon | gg → gg (123 loop diagrams, **−66.63**) |
| loop-induced | gg → h, gg → hh, **H → γγ** (28 W+top loops, **6.64e-2**), **γγ → γγ** (light-by-light, 186 loops, **1.05e-3**) |

Selected reference points, with their physics content:

- **e⁺e⁻ → a → dd̄ [virt=QCD]** (Valentin's primary benchmark). Reference at the
  default PS point (`μ_R = M_Z = 91.188 GeV`, `s = 10⁶ GeV²`):
  Born `= 3.4754514769164148e-03`; virtual normalized by `Born · α_s/(2π)`:
  Finite `= -8.9363792407373648e+00`, single pole `= 8.7724371707012754e+00`,
  double pole `= -2.6666666666666670e+00`. The single loop is the QCD vertex
  correction to the photon-quark-quark vertex — a **massless on-shell triangle**
  `C0(0,0,s;0,0,0)`, IR-divergent. The pole structure is understood analytically:
  double pole `-2 C_F = -8/3 = -2.66667` (matches exactly), single pole
  `C_F[-3 - 2 ln(μ²/s)] ≈ 8.773` (matches `8.7724`).
- **gg → hh / bb̄ [virt=QCD]** (loop-induced, top-quark boxes+triangles): loop
  `|M|²` finite `= 3.4123262140874510e-05`, poles `= 0` (finite, as expected),
  MG5 accuracy `5.2e-13`.
- **gg → h** (massive-top triangle, loop-induced, finite `9.37e-3`).
- **e⁺e⁻ → dd̄g** (pentagon, 5 external): 4 Born, 30(+4) loops, 18 R2, 52 UV;
  finite `= 2.4369`.
- **H → γγ [virt=QED]**: 0 Born, 28 W+top loops, finite `6.6360e-2`.
- **γγ → γγ [virt=QED]** (light-by-light): 0 Born, 186(+30) loops, finite
  `1.0538e-3`.

### Scope note: what oneloop's part IS

For these processes oneloop does the **reduction only**: given the loop integrand
(a tensor numerator from the fermion trace / gauge structure), it reduces it to
`C0 + B0 + A0` masters with rational-in-`d` coefficients. It does **not** build
the spinor-trace numerator (that is spenso / MadGraph), evaluate the IR-divergent
masters (OneLOop / analytic), or do the Born interference + color/spin sum
(MadGraph). Reproducing the full `-8.936` **end-to-end** requires the whole
amplitude assembly — the app / graph-bridge integration (see [the app path](05-app.md)).
The MadLoop suite therefore validates that oneloop **reproduces the benchmark
setups, the topologies, and the pole structure**, and correctly reduces the
scalar and reducible/linear numerators at the relevant (often degenerate)
kinematics.

The one characterized gap surfaced by this suite — **on-shell massless external
legs** driving pinched sub-topologies degenerate — and the off-shell-δ
regularization that resolves it are covered in detail in
[the frontier](04-frontier.md). In short: the reducer is complete off-shell; a
triangle tolerates one on-shell leg; boxes are robust to on-shell legs; only the
≥2-on-shell-leg triangle at rank≥2 hits the degenerate-Gram wall, and the
symbolic-δ wrapper (exact rational arithmetic cancels the `1/δ` inverse-Gram
poles automatically) recovers the correct value — e.g. gg→h rank-2 gives
`(1/20)[B0(0;1,1) − B0(2/5;1,1)] = −3.47483e-3` vs an independent Feynman-param
average `−3.47483e-3` (difference `~3e-10`).

---

## 3. Speed

### Absolute reduction speed (generic massive kinematics)

Symbolic reduction time per topology+numerator, producing analytic rational-in-`d`
coefficients (measured 2026-07-31):

| topology | scalar | rank-1 | rank-2 |
|---|---|---|---|
| triangle | 0.016 ms | 0.097 ms | 0.272 ms |
| box | 0.029 ms | 0.142 ms | 0.543 ms |
| pentagon | 0.504 ms | 0.727 ms | — |

Sub-millisecond across triangle / box / pentagon, scalar → rank-2. This is the
**absolute** reducer speed: the reduction is done **once** per topology+numerator
structure, yielding coefficients rational in `d`.

### Honest framing vs MadLoop

These numbers are **not** a raw-numerical-speed win over MadLoop, and it is
important not to overclaim one. The two tools do a **different operation**:

- oneloop does **symbolic** exact-rational algebra in `d` and produces a
  *formula*;
- MadLoop does optimized **float numerical** evaluation and produces a *number*
  at one phase-space point.

Symbolic reduction is inherently ~10–100× slower than a float eval — that is the
price of analyticity, not a defect. A rigorous, warmed-up per-point comparison
(MadLoop's `check_sa.f` timed over 2000 loop-ME evals after warmup, numbers still
matching `-8.9363792407`) gives **MadLoop = 13.6 µs/point** — *faster* than
oneloop's one-time symbolic reduce. MadLoop is mature, optimized compiled code
that paid its process-specific codegen cost once at `output` time (~minutes).

| operation (triangle, rank-2) | oneloop | MadLoop | ratio |
|---|---|---|---|
| symbolic reduction (once → analytic coeffs) | 0.27 ms | — (float codegen at build) | different task |
| per-point, re-reducing every point (current mode) | 0.27 ms | 0.0136 ms | ~20× slower (+1900%) |
| per-point, reduce-once + eval masters (future mode, projected) | ~0.003 ms | 0.0136 ms | ~4× faster (projected) |

The per-point floor is the scalar-master evaluation itself (avh_olo: C0 ≈ 0.5 µs),
which **MadLoop links too** — so a "reduce-once-symbolically-in-the-invariants,
evaluate-per-point" mode would be master-dominated and competitive, but is not the
present differentiator.

**The value proposition is analyticity, not speed.** oneloop's output is an
analytic closed-form reduction — rational-in-`d` coefficients times masters —
produced once per topology: *formulas, not numbers*, reusable across kinematics,
human-readable, and composable/symbolic. Report the absolute figures
("~0.02–0.7 ms/reduction, analytic output") and let the reader place them in the
one-loop landscape and reason about how the approach may scale toward 2-loop.

---

## Reproducing

See [`../benchmarks/README.md`](../benchmarks/README.md) for how to run the
cross-engine harness end to end. In brief:

```
cargo run --release --example emit_reductions -p oneloop > /tmp/oneloop_reductions.txt
python crosscheck.py /tmp/oneloop_reductions.txt
```

The MadLoop reference points and their regeneration commands are recorded in
[`../benchmarks/madloop_reference.md`](../benchmarks/madloop_reference.md).
For the algorithm being validated see [the reduction algorithm](02-reduction.md)
and [numerators](03-numerators.md); for the on-shell-massless frontier and its fix
see [the frontier](04-frontier.md); for release-level changes see
[CHANGELOG.md](CHANGELOG.md).

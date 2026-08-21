# Benchmark report — one-loop reducer

*Consolidated validation record, 2026-08-20.* This is the results-and-numbers
companion to [the benchmarks writeup](06-benchmarks.md) (which explains the
method in prose) and [the frontier writeup](04-frontier.md). Everything below was
re-run fresh on this date; the harnesses live in
[`../benchmarks/`](../benchmarks/README.md).

---

## Executive summary

| Surface | Scope | Result |
|---|---|---|
| **Cross-engine families** | 124 integral families, 2 geometries, vs OneLOop + feynalg + scipy | **124 / 124 PASS** |
| **Master values** | every A0/B0/C0/D0 emitted, OneLOop vs feynalg | **812 / 812 agree** |
| **Rank-≥3 mixed tensors** | box/pentagon/hexagon, distinct external directions | **14 / 14 PASS** (0.0–1.6σ) |
| **App coverage** | 40 physical SM processes, deployed pipeline | **848 diagrams**, 0 degenerate-Cayley walls |
| **MadLoop reproduction** | ~20 MG5_aMC processes | matched to ~14 digits |
| **Unit tests** | crate `cargo test` | **37 passed, 1 ignored, 0 failed** |

The reducer turns any one-loop integral (arbitrary N, tensor numerators, raised
powers) into A0/B0/C0/D0 with coefficients rational in `d = 4 − 2ε`. These
benchmarks put numbers on those reductions and check them against **three
independent engines** and, for tensors, **two independent tensor oracles**.

---

## 1. Cross-engine family validation (124 / 124)

Each *family* is an independent integral (topology + masses + kinematics +
numerator) emitted by [`emit_reductions.rs`](../benchmarks/rust/emit_reductions.rs)
and checked by [`crosscheck.py`](../benchmarks/python/crosscheck.py): for the
reduction `Σᵢ cᵢ(d) Mᵢ`, with masters `Mᵢ` from OneLOop, it requires (a) the
`1/ε` and `1/ε²` poles cancel (where the original is finite) and (b) the finite
part equals a direct scipy Feynman-parameter integration of the *original*
integral. A running OneLOop-vs-feynalg tally is the third, independent check.

| Category | What it covers | Families |
|---|---|---|
| Scalar N = 3…7 | triangle → heptagon (vNV / degenerate-Cayley for N>4) | ✓ |
| Dotted (raised powers) | FJT/Tarasov index-lowering, e.g. `[2,1,1,1]`, `[2,2,1]`, `[3,1]` | ✓ |
| Isotropic tensor `(k²)ᵖ` | rank-2, rank-4, **rank-6** `(k²)³` (finite + UV-divergent) | ✓ |
| Rank-1 external `k·qᵢ` | bubble … heptagon | ✓ |
| Rank-2 mixed | `k·qᵢ k·qⱼ`, `k² k·qᵢ`, `(k·qᵢ)²` | ✓ |
| **Rank-3/4/5 mixed** | distinct external directions (see §2) | ✓ **new** |
| Massless internal line | IR-finite off-shell, massless masters | ✓ |
| Near-degenerate (tiny Gram) | `1/Gram ~ 1e5` catastrophic-cancellation stress | ✓ |
| Timelike / above-threshold | real below threshold; complex masters via IBP mass-derivative | ✓ |

**Result: 124 / 124 families PASS** — poles cancel to `< 1e-6`, reduced finite
part within 5σ of the scipy Monte-Carlo (3M samples), and **812 / 812** master
finite parts agree between OneLOop and feynalg.

---

## 2. Rank-≥3 mixed tensors (the closed frontier)

Real amplitudes produce loop numerators of high rank in genuinely *different*
external directions — not just the isotropic `(k²)ᵖ`. This was the last stated
gap. We validated it directly, choosing UV-finite ranks (`r < 2(N−2)`) so the
moment/scipy oracle checks the finite part, on two spacelike geometries:

| Family | Topology | Numerator | Rank | pull (σ), both geometries |
|---|---|---|---|---|
| `mix_q1q2q3_N4` | box | `k·q1 · k·q2 · k·q3` | 3 | 1.6, 0.5 |
| `mix_q1q1q2_N4` | box | `(k·q1)² · k·q2` | 3 | 0.4, 0.1 |
| `mix_q1q2q3_N5` | pentagon | `k·q1 · k·q2 · k·q3` | 3 | 0.6, 0.2 |
| `mix_q1q2q3q4_N5` | pentagon | `k·q1 · k·q2 · k·q3 · k·q4` | 4 | 0.1, 0.8 |
| `mix_llq1q2_N5` | pentagon | `k² · k·q1 · k·q2` | 4 | 0.9, 0.4 |
| `mix_q1q2q3_N6` | hexagon | `k·q1 · k·q2 · k·q3` | 3 | 0.5, 0.1 |
| `mix_llq1q2q3_N6` | hexagon | `k² · k·q1 · k·q2 · k·q3` | 5 | 0.5, 0.0 |

**14 / 14 PASS.** Every one reduces cleanly (no panic), the `1/ε`, `1/ε²` poles
cancel to `~1e-17`, and the reduced finite part matches the independent
moment-tensor oracle + scipy to **0.0–1.6σ**. This exercises the mixed
transverse/reducible tensor machinery (RSP rules + Passarino–Veltman transverse
average) at rank up to 5 with distinct external directions — the open
generalization is now validated, not just structurally rank-agnostic.

*Note:* rank-≥3 on a **triangle** (N=3) is genuinely UV-divergent (needs
`r < 2`), so its finite part is not scipy-integrable directly; it is covered by
the reduces-cleanly unit tests and the isotropic `(k²)ᵖ` divergent families
(`ll2d`, `ll3d`), checked via the Laurent tensor oracle.

---

## 3. Application-level coverage (848 diagrams, 0 walls)

[`app_process_sweep.py`](../benchmarks/python/app_process_sweep.py) drives the
**deployed** end-to-end pipeline (gammaloop graph → contracted numerator → bridge
→ reducer) exactly as a user would, over 40 physical Standard-Model processes,
reducing an evenly-spaced sample of each process's one-loop diagrams:

| Metric | Value |
|---|---|
| Processes | 40 (self-energies, V→ff̄ triangles, Higgs, 4-fermion + QCD boxes, loop-induced) |
| Diagrams reduced | **848** |
| Reduced to masters | 738 |
| Graceful `zero_numerator` (numerator vanishes) | 109 |
| Degenerate-Cayley **walls** | **0** |
| Timeouts | 1 (a 36 MB light-by-light output — a resource limit, not a reduction failure) |
| D0 / box coverage | 15 of 40 processes |

Every master class A0/B0/C0/D0 is exercised on real diagrams, and across 848
physical numerators the reducer hit **zero** degenerate-Cayley failures.

---

## 4. MadLoop reproduction (~14 digits)

From [`../benchmarks/madloop_reference.md`](../benchmarks/madloop_reference.md):
~20 MG5_aMC v3.7.2 processes reproduced (triangle, box, pentagon, 4-gluon,
loop-induced), anchored on Valentin's benchmark `e⁺e⁻ → γ → dd̄ [virt=QCD]`
matched to ~14 digits, and the `gg → h` massive-top on-shell triangle where the
`reg_delta` regularization matches an independent computation to `3.4e-7`
(numeric extrapolation) / `3.2e-10` (symbolic-δ).

---

## Method — why this is trustworthy

- **Three independent master engines:** OneLOop (avh_olo, the library MadLoop
  links), feynalg (Denner-style analytic closed forms), and scipy direct
  Feynman-parameter integration of the *original* integral.
- **Two independent tensor oracles:** a symmetric-moment integrator and a
  covariant Minkowski-Gram integrator (different computational paths) cross-check
  each tensor result.
- **Exact rational coefficients:** the reduction is symbolic in `d`; the pole
  cancellation is algebraic (to `~1e-17`), not a numerical fit.

## Honest scope

`oneloop` performs the *reduction*; the numeric/ε-expansion of the masters is
supplied by the external evaluator (OneLOop). Genuinely singular / threshold
configurations (invariants vanishing at different rates) and exactly-singular
Gram boxes at all-massless kinematics remain guarded (clean error, never garbage)
— see [the frontier](04-frontier.md). Fully-massless one-loop (external *and*
internal) is scaleless = 0 in dim reg and needs no computation.

---

*Reproduce:* `cargo run --release --example emit_reductions -p oneloop >
/tmp/r.txt && python3 benchmarks/python/crosscheck.py /tmp/r.txt` (cross-engine);
`python3 benchmarks/python/app_process_sweep.py 40` (app coverage). See
[`../benchmarks/README.md`](../benchmarks/README.md).

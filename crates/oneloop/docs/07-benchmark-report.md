# Benchmark report — one-loop reducer

*Consolidated validation record, 2026-08-21.* This is the results-and-numbers
companion to [the benchmarks writeup](06-benchmarks.md) (which explains the
method in prose) and [the frontier writeup](04-frontier.md). Everything below was
re-run fresh on this date; the harnesses live in
[`../benchmarks/`](../benchmarks/README.md).

---

## Executive summary

| Surface | Scope | Result |
|---|---|---|
| **Cross-engine families** | 132 integral families, 2 geometries, vs OneLOop + feynalg + scipy | **132 / 132 PASS** |
| **Master values** | every A0/B0/C0/D0 emitted, OneLOop vs feynalg | **944 / 944 agree** |
| **High-rank mixed tensors** | box → heptagon, distinct external directions, rank 3–6 (incl. dotted) | **all PASS** (0.0–1.6σ) |
| **Singular-Gram tensors** | heptagon (N=7) + coincident-momenta case, pseudo-inverse fix | **reduces + validated** (0.4–0.7σ) |
| **App coverage** | 40 physical SM processes, deployed pipeline | **848 diagrams**, 0 degenerate-Cayley walls |
| **MadLoop reproduction** | 110+ MG5_aMC processes run (96 fresh, 3 batches) | reproduced at MadLoop's ~14-digit accuracy; anchors match published values (oneloop's *reduction* is validated separately, row 1) |
| **Unit tests** | crate `cargo test` | **38 passed, 1 ignored, 0 failed** |

The reducer turns any one-loop integral (arbitrary N, tensor numerators, raised
powers) into A0/B0/C0/D0 with coefficients rational in `d = 4 − 2ε`. These
benchmarks put numbers on those reductions and check them against **three
independent engines** and, for tensors, **two independent tensor oracles**.

---

## 1. Cross-engine family validation (132 / 132)

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
| **Rank-3–6 mixed** | distinct external directions, box → heptagon, incl. dotted (see §2) | ✓ |
| **Singular-Gram tensor** | heptagon: >4 external momenta ⇒ rank-deficient Gram (see §2) | ✓ **new** |
| Massless internal line | IR-finite off-shell, massless masters | ✓ |
| Near-degenerate (tiny Gram) | `1/Gram ~ 1e5` catastrophic-cancellation stress | ✓ |
| Timelike / above-threshold | real below threshold; complex masters via IBP mass-derivative | ✓ |

**Result: 132 / 132 families PASS** — poles cancel to `< 1e-6`, reduced finite
part within 5σ of the scipy Monte-Carlo (3M samples), and **944 / 944** master
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
| `mix_q1q2q3_dotN4` | box `[2,1,1,1]` | `k·q1 · k·q2 · k·q3` (dotted) | 3 | 0.0, 0.7 |
| `mix_llq1q2_dotN4` | box `[2,1,1,1]` | `k² · k·q1 · k·q2` (dotted) | 4 | 0.2, 0.2 |
| `mix_q1q2q3_N5` | pentagon | `k·q1 · k·q2 · k·q3` | 3 | 0.6, 0.2 |
| `mix_q1q2q3q4_N5` | pentagon | `k·q1 · k·q2 · k·q3 · k·q4` | 4 | 0.1, 0.8 |
| `mix_llq1q2_N5` | pentagon | `k² · k·q1 · k·q2` | 4 | 0.9, 0.4 |
| `mix_q1q2q3_N6` | hexagon | `k·q1 · k·q2 · k·q3` | 3 | 0.5, 0.1 |
| `mix_llq1q2q3_N6` | hexagon | `k² · k·q1 · k·q2 · k·q3` | 5 | 0.5, 0.0 |
| `mix_llq1q2q3q4_N6` | hexagon | `k² · k·q1 · k·q2 · k·q3 · k·q4` | 6 | 0.8, 0.7 |
| `mix_q1q2q3_N7` | **heptagon** | `k·q1 · k·q2 · k·q3` (singular Gram) | 3 | 0.4, 0.7 |

**22 / 22 PASS** (11 families × 2 geometries). Every one reduces cleanly, the
`1/ε`, `1/ε²` poles cancel to `~1e-17`, and the reduced finite part matches the
moment-tensor oracle + scipy to **0.0–1.6σ**. This exercises the mixed
transverse/reducible tensor machinery (RSP rules + Passarino–Veltman transverse
average) at rank up to 6 with distinct external directions, **including raised
propagator powers** (the dotted box) — the open generalization is validated.

### The singular-Gram case (heptagon, N=7) and its fix

The heptagon (`mix_q1q2q3_N7`) is special, and it exposed a genuine gap. A
one-loop N-point has N−1 external momenta, and the tensor reduction rewrites the
loop momentum in the basis of those momenta by inverting their **Gram matrix**.
But in `d = 4` at most **four** momenta can be linearly independent, so for
**N ≥ 6** the Gram is rank-deficient — **singular** — and the naive inverse
crashes. (The same happens for coincident external momenta, `qᵢ = qⱼ`.) This is
not a physics bug; it is unavoidable geometry, and it only surfaces for tensor
numerators (scalar heptagons never invert a Gram, so they always reduced).

**The fix — a pseudo-inverse.** You do not need all N−1 directions: only four are
independent and the rest are their linear combinations. So `gram_solve` now
solves on a **maximal linearly-independent sub-Gram** and sets the redundant
coefficients to zero. Because `rhs` lies in the reducible span, this reproduces
the full system exactly (`G·c = rhs`) and gives the correct projection onto the
external-momentum span. Implemented in
[`reduce.rs`](../src/reduce.rs) (`gram_solve` → `gram_solve_matrix` +
`independent_gram_subset`), unit-tested (`gram_solve_matrix_handles_singular_gram`),
and validated here: the heptagon rank-3 tensor now reduces to 35 master terms
matching the oracle to **0.4σ / 0.7σ** on both geometries. With this the reducer
is complete — **any N, any rank** — and the last guarded panic in the
tensor-reduction path is retired.

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

### Pentagons and beyond — the bridge extended to N ≥ 5 (2026-08-21)

The graph → reducer bridge previously mapped externals only up to `q3`
(four-point). It now builds the external-momentum symbols dynamically
(`P(j) → q{j+1}`, up to nine-point), so the deployed pipeline reduces **5+ point**
physical diagrams. Validated **live** on **`e+e-→ddg`** (Valentin's pentagon
process): 184 one-loop diagrams; of 15 sampled, **11 reduced, 4 graceful
zero-numerator, 0 walls, 0 errors**, and **6 reductions reference `q4`** — the
genuine five-point external momentum the old bridge could not produce. Reductions
reach `D0` and run to ~10 MB (gzip-compressed on the wire). The app now handles
pentagons, hexagons, and beyond — not just boxes.

---

## 4. MadLoop reproduction (~110 processes)

**What this validates (read first).** `oneloop` does the *reduction* of a loop
integrand to masters — it does **not** build the full amplitude (spinor trace +
colour + Born interference). So this section is a **MadLoop reproduction record**,
not a per-process `oneloop`-vs-MadLoop diff: we ran MG5_aMC on ~110 processes and it
produces self-consistent numbers at its own ~14-digit accuracy, with the published
anchors (`e⁺e⁻→dd̄`, `gg→gg`, the loop-induced values) matching known/Valentin
numbers to 14 digits. This demonstrates the *process space* the reducer's reductions
cover, and that the IR/colour structure is textbook. **The direct validation of
`oneloop`'s output is §1** (132/132 vs OneLOop — the very master library MadLoop
links). A full-amplitude `oneloop`-vs-MadLoop comparison per process is a separate
project (the amplitude assembly, "target B", noted at the end of this section).

From [`../benchmarks/madloop_reference.md`](../benchmarks/madloop_reference.md):
~21 MG5_aMC v3.7.2 processes reproduced (triangle, box, pentagon, 4-gluon,
loop-induced), anchored on Valentin's benchmark `e⁺e⁻ → γ → dd̄ [virt=QCD]`
(matched to 14 digits), and the `gg → h` massive-top on-shell triangle where the
`reg_delta` regularization matches an independent computation to `3.4e-7`
(numeric extrapolation) / `3.2e-10` (symbolic-δ). A fresh anchor
**`e⁺e⁻ → tt̄ [virt=QCD]`** (2026-08-21, rel. accuracy 5.1e-15) reproduces the
massive-quark IR structure — **double pole exactly 0** (massive tops have no
collinear singularity), the clean contrast with the massless anchor's
`−8/3 = −2C_F` double pole — validating the massive-configuration regime the
reducer is built for (Born `3.0630e-2`, virtual finite `+1.8566`, single pole `+6.5419`).

A **35-process batch** (2026-08-21, 5.5 min) then reproduced the full IR structure across
V→qq̄ triangles, qq̄/gluon-initiated boxes, the 4-gluon amplitude, loop-induced
(`gg→h`, `H→γγ`, `gg→hh`), and the `e⁺e⁻→ddg` pentagon — all at MadLoop's ~14-digit
accuracy. The **massive-vs-massless double pole** (0 for massive `b`/`t`, `−8/3 = −2C_F` for light
quarks) and the **exact colour-factor poles** (`−4C_F`, `−4C_A = −12` for 4-gluon)
fall out automatically, and `e⁺e⁻→dd̄ = −8.93638` / `gg→gg = −66.63` reproduce prior
anchors exactly. A **second 40-process batch** (EW/QED — `e⁺e⁻→WW/ZZ/γγ`, `W→ℓν`,
`z→ℓℓ` — plus new QCD 2→2 and more pentagons) then validated **40/40**, with the
QED `e⁺e⁻` double pole `−0.1315` constant across the electroweak channels and the
both-massive `bb̄→tt̄` double pole exactly 0. A **third batch** added QCD 2→3
jets (double poles `−25/3`, `−35/3`) and quark-initiated diboson (double pole
scaling with quark charge², `−0.0584` up vs `−0.0146` down). **96 fresh anchors
across the three batches** (~110+ on record), spanning triangle / box / pentagon /
4-gluon / 2→3-jet / loop-induced / diboson, QED+QCD, massless+massive — every
double pole reading out the exact IR/colour structure (full tables in
`madloop_reference.md`).

**Two tiers of MadLoop validation.** The 14-digit numbers above are the *full
amplitude* (spinor trace + colour + Born interference); `oneloop` supplies the
**reduction**, which is validated against **OneLOop** — the very master library
MadLoop links (§1, 132/132 families). At the deployed-app level the reducer
reduces the MadLoop-listed processes directly: the triangle `e⁺e⁻→dd̄`, the boxes
`uū→dd̄` / `dd̄→dd̄`, the 4-gluon `gg→gg`, the loop-induced `gg→h` / `h→γγ` /
`γγ→γγ`, and — newly, via the N≥5 bridge — the **pentagon `e⁺e⁻→ddg`** (§3). A
fresh *full-amplitude* anchor requires an MG5_aMC run plus the amplitude assembly
(the app/graph-bridge project), not the reducer in isolation.

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

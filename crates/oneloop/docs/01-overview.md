# oneloop — Overview

`oneloop` is a symbolic one-loop IBP (integration-by-parts) reducer. Given any
one-loop Feynman integral with an arbitrary polynomial numerator, it reduces the
integral to a linear combination of the four standard scalar master integrals

| master | topology | emitted head |
|--------|----------|--------------|
| A0 | tadpole  | `A0(m²)` |
| B0 | bubble   | `B0(p², m1², m2²)` |
| C0 | triangle | `C0(p1², p2², p12², m1², m2², m3²)` |
| D0 | box      | `D0(p1², p2², p3², p4², s, t, m1²..m4²)` |

with coefficients that are **rational functions of the dimension `d = 4 − 2ε`**.
The argument order of the emitted heads follows the AVH / OneLOopBridge
convention (see `src/masters.rs`).

The reduction uses **closed-form, per-topology IBP recursions** implemented in
[Symbolica](https://symbolica.io). There is **no Laporta engine**: the general
Laporta linear-system approach is deliberately scoped to a future ≥2-loop
effort, not built here.

## Goal and scope

The one sentence: **reduce, do not evaluate.**

`oneloop` is the *REDUCE* step of a two-step pipeline:

1. **REDUCE** (this crate): integral + numerator → `Σ cᵢ · Mᵢ`, where each `cᵢ`
   is a symbolic rational-in-`d` coefficient and each `Mᵢ` is one of A0/B0/C0/D0
   carried as an *opaque* function atom. No analytic form for the masters is
   baked into the core — by design.
2. **EVALUATE** (delegated): put numbers on the master symbols. This is handed
   off to **OneLOopBridge** (numeric, the `avh_olo` library) and **feynalg**
   (analytic), driven from the `benchmarks/` harnesses — not from the crate
   itself.

Concretely, `reduce(&IntegralFamily) -> Reduction` returns a `Reduction` whose
`terms` field is a `Vec<(Atom, MasterIntegral)>` of coefficient/master pairs (`src/reduce.rs`), and
`amplitude(&IntegralFamily) -> Atom` folds those into a single symbolic
expression `Σ cᵢ · symbol(Mᵢ)` (`src/amplitude.rs`). Nothing in the core assigns
a numeric or analytic value to A0/B0/C0/D0.

### Why "reduce, not evaluate" is the whole value proposition

The differentiator is **analytic closed-form output**, not raw numerical speed.
The reducer produces formulas (rational-in-`d` coefficients × master symbols)
that are reusable across kinematics, human-readable, and composable — as opposed
to a single float at a single phase-space point. Benchmarks bear this out:
per-point *numerical* evaluation in a mature tool like MadLoop is faster than
`oneloop`'s *symbolic* reduction (MadLoop's warmed-up per-point loop-ME eval is
13.6 µs; a `oneloop` triangle rank-2 symbolic reduction is 0.27 ms). We do **not**
claim a per-point speed advantage over MadLoop. See
[the benchmarks](06-benchmarks.md) for the honest speed story.

## Dispatch at a glance

`reduce()` dispatches on the number of propagators `N` (`src/reduce.rs`):

| N | reduction |
|---|-----------|
| 1 | tadpole → A0 (`tadpole_coefficient` / `tadpole_numerator`) |
| 2 | bubble → B0 + tadpoles (`reduce_cayley` / `bubble_numerator`) |
| 3 | triangle → C0 + sub-topologies (`reduce_cayley` / `triangle_numerator`) |
| 4 | box → D0 + sub-topologies (`reduce_cayley` / `box_numerator`) |
| N>4 | van Neerven–Vermaseren bordered-Cayley reduction to boxes (`reduce_cayley` / `ngon_numerator`) |

Scalar integrals of any `N` go through the modified-Cayley matrix
`Yᵢⱼ = mᵢ² + mⱼ² − (rᵢ − rⱼ)²` and its bordered determinant; raised propagator
powers use FJT/Tarasov index-lowering (`reduce_cayley`); tensor numerators go
through a Passarino–Veltman-style transverse projection and RSP substitution
(`reduce_num`). The mechanics of all of this live in
[the reduction algorithm](02-reduction.md) and [numerators](03-numerators.md).

## Current status

Grounded in the validation record (`benchmarks/madloop_reference.md`,
`benchmarks/MONDAY_AGENDA.md`).

### What works

- **Scalar integrals, any `N`** (N=1 up, verified scalar N=3–7 exact vs
  OneLOop/feynalg on generic kinematics).
- **Tensor numerators** up to rank-6 on generic massive kinematics; rank-1 and
  rank-2 are validated broadly, including on triangle and box.
- **Raised propagator powers** (e.g. exponents `[2,1,1]`) via FJT/Tarasov
  index-lowering.
- **N>4 scalar** via the van Neerven–Vermaseren bordered-Cayley reduction to
  boxes.
- **On-shell massless external legs** via off-shell `reg_delta`
  regularization: when any invariant `is_zero()`, `reduce()` replaces the
  vanishing invariants with a symbol `oneloop::reg_delta`, calls the core
  reducer (whose exact rational arithmetic makes the `1/δ` inverse-Gram poles
  explicit and cancels them in the sum `Σ cᵢ Mᵢ`), then substitutes `δ → 0`.
  This is a thin wrapper — **zero changes to the core reducer**
  (`reduce_regularized` in `src/reduce.rs`). Validated: gg→h massive-top
  triangle, rank-2, on-shell massless legs, matches an independent
  Feynman-parameter computation to 3.4e-7 (numeric) / 3e-10 (symbolic-δ).

### Validated against

- **MadLoop** (MG5_aMC v3.7.2): **~110 processes reproduced** (96 fresh, 3 batches)
  at MadLoop's ~14-digit accuracy; anchors match published values to 14 digits —
  including e+e⁻→γ→dd̄ `[virt=QCD]` (triangle, finite part −8.936),
  gg→hh (loop-induced box), gg→h (massive-top triangle), e+e⁻→dd̄g (pentagon),
  Drell-Yan, W/Z production, 4-gluon gg→gg (123 loop diagrams, −66.63), and
  loop-induced H→γγ (28 W+top loops, 6.64e-2) and γγ→γγ (186 loops, 1.05e-3).
- **OneLOopBridge + feynalg + scipy** cross-engine harness (under
  `benchmarks/`): scalar N=3–7 exact on generic kinematics; tensor to rank-6;
  UV-divergent, massless-line, timelike/threshold/near-degenerate points.

### Known gaps / frontier

The reducer is a fundamentally **massive-configuration** method (mass-derivative
/ Tarasov recursions that divide by Gram and Cayley determinants). It is
**complete off-shell**; the frontier is where those determinants genuinely
vanish:

- **Massless internal propagators** (m²=0, IR-divergent regime) — the
  `reg_delta` wrapper covers many of these (e.g. it cures the previously-panicking
  e+e⁻→dd̄ `dot(k, q₂)` case), but genuinely singular / threshold configurations
  where different invariants must vanish at *different rates* are not covered by
  the single-shared-`δ` limit.
- **Exactly-singular Gram** configurations that survive regularization are
  **guarded** — the reducer returns a clean error (`OneLoopError`,
  `src/error.rs`) or panics with a specific message rather than emitting garbage.

The rigorous general solution to the residual cases is the Denner–Dittmaier /
Tarasov dimension-shift `d → d+2` machinery. This is characterized in detail in
[the frontier](04-frontier.md).

## How it fits gammaloop

`oneloop` is a workspace crate inside `gammaloop-alphaloop`. It consumes a
one-loop graph produced by gammaloop and turns it into master integrals:

```
gammaloop Graph ──▶ src/bridge.rs ──▶ IntegralFamily ──▶ reduce() ──▶ Σ cᵢ · Mᵢ
                    (numerator_to_dot_form,               (masters,
                     family_from_gammaloop)                rational in d)
```

`src/bridge.rs` translates a gammaloop one-loop graph — its contracted numerator
atom and its per-edge loop-momentum-basis representation + masses — into an
`IntegralFamily` the reducer can consume (`numerator_to_dot_form`,
`external_offset_from_lmb_rep`, `family_from_gammaloop`). Symbolica's license is
activated through `gammalooprs`; the crate's *tests* use the `#[cfg(test)]`
helper `ensure_symbolica_license()` in `src/lib.rs`. The end-to-end app
integration is covered in [the app](05-app.md).

## File map (`src/`)

| file | responsibility |
|------|----------------|
| `lib.rs` | crate root; module wiring; the test-only `ensure_symbolica_license()` helper (`#[cfg(test)]`; activates the Symbolica license via `gammalooprs`) |
| `error.rs` | `OneLoopError` — `UnsupportedLoopOrder`, `ExtractionFailed`, `Symbolica` |
| `family.rs` | `IntegralFamily`, `Propagator`, `Kinematics`, `Integral`, `Isp` — the input data model |
| `symbols.rs` | registered Symbolica symbols behind the `S` singleton: `d`, `k`, `q1..q3`, `psq`, `dot` (symmetric + linear), master heads `A0/B0/C0/D0` |
| `masters.rs` | `MasterIntegral` enum + `MasterBasis::symbol()` emitter (opaque master atoms, AVH arg order) |
| `amplitude.rs` | `amplitude(&IntegralFamily) -> Atom`: folds `Σ cᵢ · symbol(Mᵢ)` |
| `reduce.rs` | `reduce()` dispatcher, `reduce_regularized` wrapper, and every per-topology reducer (Cayley, FJT/Tarasov, PV projection) |
| `bridge.rs` | gammaloop → oneloop translation (`numerator_to_dot_form`, `family_from_gammaloop`) |

## Read next

- [02-reduction.md](02-reduction.md) — the reduction algorithm: modified-Cayley
  matrix, bordered determinants, FJT/Tarasov index-lowering, N>4 recursion.
- [03-numerators.md](03-numerators.md) — how arbitrary polynomial numerators in
  `dot(k, k)` and `dot(k, qᵢ)` are handled (PV transverse projection, RSP
  substitution, pinching).
- [04-frontier.md](04-frontier.md) — the on-shell / massless / degenerate-Gram
  frontier, the `reg_delta` off-shell regularization, and where the rigorous fix
  lives.
- [05-app.md](05-app.md) — the gammaloop graph bridge and end-to-end app path.
- [06-benchmarks.md](06-benchmarks.md) — validation record and the honest speed
  story vs MadLoop.
- [../benchmarks/README.md](../benchmarks/README.md) — the cross-check harness.
- [CHANGELOG.md](CHANGELOG.md) — change history.

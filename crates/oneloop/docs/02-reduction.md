# The reduction algorithm

This document describes how `oneloop` reduces a one-loop integral family to the
four scalar master integrals. It covers the master-integral data model and the
`symbol()` emitter, the per-topology closed-form IBP recursions (`N = 0…4`), the
`N > 4` reduction to boxes via van Neerven–Vermaseren, and the `reduce()`
dispatch table. Source files: [`src/reduce.rs`](../src/reduce.rs),
[`src/masters.rs`](../src/masters.rs), [`src/family.rs`](../src/family.rs).

For the crate scope and where this step sits, see
[the overview](01-overview.md); for how arbitrary polynomial numerators are
handled, see [numerators](03-numerators.md); for the kinematic regimes that are
still guarded, see [the frontier](04-frontier.md).

---

## The four scalar masters

Every one-loop integral in `d = 4 − 2ε` reduces to a linear combination of four
scalar master integrals, with coefficients rational in `d`:

| Master | Points | Head atom | Emitter arguments (AVH / OneLOopBridge order) |
|--------|--------|-----------|-----------------------------------------------|
| A0 (tadpole)  | 1 | `oneloop::A0` | `A0(m²)` |
| B0 (bubble)   | 2 | `oneloop::B0` | `B0(p², m1², m2²)` |
| C0 (triangle) | 3 | `oneloop::C0` | `C0(p1², p2², p12², m1², m2², m3²)` |
| D0 (box)      | 4 | `oneloop::D0` | `D0(p1², p2², p3², p4², s, t, m1², m2², m3², m4²)` |

The core reducer produces **opaque** master atoms only — it never inlines an
analytic ε-expansion. Putting numbers on A0/B0/C0/D0 is delegated to an
evaluator (OneLOopBridge numerically, feynalg analytically); see
[benchmarks](06-benchmarks.md) and [the app](05-app.md).

### The `MasterIntegral` enum

Masters are carried as a typed enum whose variant fields are Symbolica `Atom`s
holding the kinematic invariants and squared masses ([`masters.rs`](../src/masters.rs)):

```rust
pub enum MasterIntegral {
    Tadpole  { m_sq: Atom },
    Bubble   { p_sq: Atom, m1_sq: Atom, m2_sq: Atom },
    Triangle { p1_sq, p2_sq, p12_sq, m1_sq, m2_sq, m3_sq: Atom },
    Box      { p1_sq, p2_sq, p3_sq, p4_sq, s, t, m1_sq, m2_sq, m3_sq, m4_sq: Atom },
}
```

### The `symbol()` emitter and argument order

The `MasterBasis` trait maps a `MasterIntegral` to its head-function atom.
`OneLoopMasters` is the concrete basis; its `symbol()` is total (infallible) and
fixes the argument order to match the AVH / OneLOopBridge convention so the
emitted atoms feed the numeric/analytic evaluators directly:

```rust
impl MasterBasis for OneLoopMasters {
    fn symbol(&self, integral: &MasterIntegral) -> Atom {
        match integral {
            Tadpole  { m_sq }                       => function!(S.a0, m_sq),
            Bubble   { p_sq, m1_sq, m2_sq }         => function!(S.b0, p_sq, m1_sq, m2_sq),
            Triangle { p1_sq, p2_sq, p12_sq, .. }   => function!(S.c0, p1_sq, p2_sq, p12_sq, m1_sq, m2_sq, m3_sq),
            Box      { p1_sq, .., m4_sq }           => function!(S.d0, p1_sq, p2_sq, p3_sq, p4_sq, s, t, m1_sq, m2_sq, m3_sq, m4_sq),
        }
    }
}
```

The head symbols `S.a0 = oneloop::A0`, `S.b0 = oneloop::B0`,
`S.c0 = oneloop::C0`, `S.d0 = oneloop::D0` are registered in
[`symbols.rs`](../src/symbols.rs). The `d`-dependence lives entirely in the
scalar coefficients (via the registered symbol `S.d = oneloop::d`), not in the
master arguments.

---

## The data model: `IntegralFamily`

The reducer consumes an `IntegralFamily` ([`family.rs`](../src/family.rs)):

```rust
pub struct IntegralFamily {
    pub propagators: Vec<Propagator>,   // { momentum: Atom, mass_sq: Atom }
    pub isps:        Vec<Isp>,          // reserved (currently unused)
    pub kinematics:  Kinematics,        // { invariants: Vec<Atom> }
    pub targets:     Vec<Integral>,     // { propagator_exponents: Vec<i32>, .. }
    pub numerator:   Atom,              // polynomial in oneloop::dot(k, ·)
}
```

- **`propagators`** — one `Propagator` per loop line; `mass_sq` is the squared
  mass `m_i²`. The propagator count `family.propagators.len()` is what the
  dispatcher branches on.
- **`kinematics.invariants`** — the `C(N, 2)` pairwise invariants `(r_i − r_j)²`
  in lexicographic `(i < j)` order (`r_i` is the external offset of line `i`).
  A regular `N`-point family therefore carries `N(N−1)/2` invariants; the
  `N > 4` branch asserts exactly this count.
- **`targets[0].propagator_exponents`** — the propagator powers `a_i` of the
  target topology (all `1` for a plain scalar integral; `> 1` for dotted
  propagators).
- **`numerator`** — a polynomial in the symmetric linear symbol
  `oneloop::dot`, built from `dot(k, k)` and `dot(k, q_i)` for the external
  momenta `q1, q2, q3` (registered in [`symbols.rs`](../src/symbols.rs)). A
  purely scalar integral has `numerator == Atom::num(1)`. Numerator handling is
  the subject of [numerators](03-numerators.md).

---

## `reduce()` dispatch

The public entry point is:

```rust
pub fn reduce(family: &IntegralFamily) -> Reduction
```

returning a `Reduction { terms: Vec<(Atom, MasterIntegral)> }` — a list of
(coefficient, master) pairs. `Reduction::simplify()` cancels each coefficient to
lowest terms.

`reduce()` is a thin wrapper. If **any** kinematic invariant `is_zero()`
(on-shell massless external legs, where modified-Cayley / Gram determinants can
vanish), it routes to `reduce_regularized()`, which replaces every zero
invariant with a symbol `oneloop::reg_delta`, calls `reduce_core()`, then
substitutes `δ → 0` in the resulting coefficients and master arguments and drops
zero terms. The `1/δ` inverse-Gram poles cancel algebraically in the sum
`Σ cᵢ Mᵢ`, so this off-shell regularization is exact. Otherwise it calls
`reduce_core()` directly.

`reduce_core()` dispatches on the propagator count:

| Props `N` | Scalar (`numerator == 1`) | With numerator |
|-----------|---------------------------|----------------|
| 0 | empty (scalar `1`) | — |
| 1 | `tadpole_coefficient` → A0 | `tadpole_numerator` |
| 2 | modified-Cayley → B0 (+ tadpoles) | `bubble_numerator` |
| 3 | modified-Cayley → C0 | `triangle_numerator` |
| 4 | modified-Cayley → D0 | `box_numerator` |
| `> 4` | `modified_cayley` + `reduce_cayley` (→ boxes) | `ngon_numerator` |

The scalar branch is uniform: for `N ≥ 2` it builds the modified Cayley matrix
and hands it to `reduce_cayley`, so all of B0/C0/D0/`N`-gon scalar reduction go
through one code path. The per-topology *numerator* routines (`bubble_numerator`,
`triangle_numerator`, `box_numerator`, `ngon_numerator`) are covered in
[numerators](03-numerators.md).

---

## Per-topology closed-form IBP recursions

### N = 0 — scalar 1

No propagators: the reduction is empty (the bare scalar `1`).

### N = 1 — tadpole → A0

A single line with power `a₁` and mass `m²` reduces to A0 with a
closed-form rational coefficient (`tadpole_coefficient`):

```
coeff(a₁, m²) = Π_{k=1}^{a₁−1}  (d − 2k) / (2k · m²)          for a₁ ≥ 1,
              = 0                                              if a₁ ≤ 0 or m² = 0.
```

So `a₁ = 1` gives coefficient `1` (the plain tadpole), and each raised power
peels off a `(d − 2k)/(2k·m²)` factor. A massless tadpole (`m² = 0`) or a
non-positive power is scaleless and returns `0`.

### N = 2 — bubble → B0 (+ tadpoles)

For a scalar bubble the reducer builds the `2×2` modified Cayley matrix from the
two masses and the single invariant `p²` and calls `reduce_cayley`. With unit
powers this emits `B0(p², m1², m2²)` directly. Raised propagator powers are
lowered by the FJT/Tarasov index-lowering step inside `reduce_cayley` (below),
producing B0 plus the two single-mass tadpoles `A0(m1²)`, `A0(m2²)`.

### N = 3 — triangle → C0

The scalar triangle builds the `3×3` modified Cayley matrix from the three
masses and the three invariants and reduces via `reduce_cayley`, emitting
`C0(p1², p2², p12², m1², m2², m3²)`. The three invariants fill the lexicographic
Cayley slots in order — `(s₀₁, s₀₂, s₁₂) = (inv(0), inv(1), inv(2))`, no remap —
while the emitted C0 arguments are the permutation
`(p1², p2², p12²) = (inv(0), inv(2), inv(1)) = (s₀₁, s₁₂, s₀₂)`.

### N = 4 — box → D0

The scalar box builds the `4×4` modified Cayley matrix from four masses and six
invariants and reduces via `reduce_cayley`, emitting the ten-argument
`D0(p1², p2², p3², p4², s, t, m1², m2², m3², m4²)`. The six physics invariants
are remapped to the lexicographic Cayley order
`(s₀₁, s₀₂, s₀₃, s₁₂, s₁₃, s₂₃) = (p1, s, p4, p2, t, p3)`.

### The shared engine: modified Cayley + `reduce_cayley`

The `N = 2…4` scalar reductions and the `N > 4` recursion all run through the
same two functions.

**Modified Cayley matrix.** `modified_cayley(masses, pairwise)` builds

```
Y_ij = m_i² + m_j² − (r_i − r_j)²        (i ≠ j),
Y_ii = 2 m_i²,
```

from the masses and the `C(N,2)` pairwise invariants in lexicographic order.

**`reduce_cayley(y, masses, exponents)`** is the recursion that turns a Cayley
matrix + propagator powers into masters:

1. `n = 0` → empty; `n = 1` → tadpole via `tadpole_coefficient`.
2. If any exponent is `0`, that line is pinched: drop its row/column and recurse
   on the sub-topology.
3. If all exponents are `1` and `n ≤ 4`, emit the leaf master directly
   (`emit_master`: `n = 2` → Bubble, `n = 3` → Triangle, `n = 4` → Box). The
   `emit_master` step reconstructs each `pᵢⱼ² = mᵢ² + mⱼ² − Yᵢⱼ` from the Cayley
   entries, in the same physics ordering used by `symbol()`.
4. If all exponents are `1` and `n > 4`, apply the van Neerven–Vermaseren step
   (below).
5. Otherwise (raised powers), lower the index with the FJT/Tarasov recursion:
   pick the line `k` with the largest power, decrement it, and expand using the
   inverse Cayley matrix. The coefficients carry the `d`-dependence through
   factors `(d − (total + aᵢ))` divided by the decremented power — this is where
   the rational-in-`d` structure enters. A vanishing `det(Y)` (degenerate
   kinematics) is handled by a separate degenerate branch (`degenerate_coeffs`)
   that solves the bordered null-space instead of inverting.

---

## N > 4 — van Neerven–Vermaseren reduction to boxes

In four dimensions a scalar `N`-point integral with `N > 4` is not independent:
it is a linear combination of `N` pinched `(N−1)`-point integrals. `oneloop`
implements this with the **bordered Cayley** construction (`high_point_coeffs`):

```
c_i = (adj YB)_{0i} / det(Y),
```

where `YB` is the Cayley matrix `Y` bordered by a row and column of `1`s, and
the `c_i` are the coefficients of the `N` sub-topologies obtained by pinching
line `i`. `reduce_cayley` then recurses on each `delete_row_col(y, i, i)` with
unit powers, so an `N`-gon cascades down through pentagons/boxes until every
term is a genuine four-point (or lower) master. Scalar high-point families thus
reduce to a sum of **boxes** (plus lower masters when a sub-Cayley is
degenerate).

When `det(Y) = 0` (exceptional kinematics), `high_point_coeffs` falls back to
`degenerate_coeffs`, which reads the reduction coefficients off the null-space of
the bordered matrix.

This is validated: `reduce()` on a generic massive **pentagon**
(`[1,2,3,4,5]` masses, generic invariants) returns exactly five boxes with the
van Neerven–Vermaseren coefficients

```
[ −1088/639,  −212/639,  −44/213,  −191/639,  −341/639 ],
```

and the last term is the box on propagators 1–4,
`D0(3, 4, 5, 7, 5, 6; 1, 2, 3, 4)` (see the `scalar_pentagon_reduces_to_five_boxes`
test in [`reduce.rs`](../src/reduce.rs)). Hexagons and heptagons recurse the same
way down to boxes.

---

## What is validated, and where the frontier is

- **Scalar `N = 1…N`** (any propagator count) reduces exactly. Generic-kinematics
  scalar `N = 3…7` is cross-checked against OneLOopBridge + feynalg; see
  [benchmarks](06-benchmarks.md).
- The reducer reproduces MadLoop (MG5_aMC v3.7.2) across ~20 processes, including
  the massless on-shell triangle from `e⁺e⁻ → γ → dd̄ [virt=QCD]` matched to ~14
  digits, and the `gg → h` massive-top triangle where the on-shell massless
  regularization matches an independent computation to `3.4e-7` (numeric
  extrapolation) / `3e-10` (symbolic-δ).
- **On-shell massless legs — handled by `reg_delta`.** Where `det(Y) → 0` from
  on-shell massless external legs, `reduce()` regularizes the vanishing invariants
  off-shell and takes the `δ → 0` limit exactly (the `1/δ` inverse-Gram poles
  cancel in `Σ cᵢ Mᵢ`). This covers the case the bare massive method cannot do —
  the **multi-on-shell-leg triangle at rank ≥ 2** (unit test
  `on_shell_massless_triangle_rank2_regularizes`; live across `a→dd̄`,
  `e⁺e⁻→dd̄`, `z→dd̄` on the deployed app). See [the frontier](04-frontier.md).
- **Remaining gaps** (guarded — clean panic at `reduce.rs`, never garbage):
  genuinely singular / **threshold** configurations where different invariants
  must vanish at *different rates* (the single-shared-δ limit does not cover these;
  they need the Denner–Dittmaier systematic expansion), and exactly-singular Gram
  boxes at all-massless kinematics. IR-divergent masters (`B0(0;0,0)`, scaleless
  `C0`) are reduced *to* but evaluated by the external master library, not
  `oneloop`. See [the frontier](04-frontier.md) for the full map.

---

See also: [../benchmarks/README.md](../benchmarks/README.md) for the validation
harness and [CHANGELOG.md](CHANGELOG.md) for the reduction-feature history.

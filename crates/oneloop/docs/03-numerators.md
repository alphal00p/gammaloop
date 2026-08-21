# Tensor & dotted numerators

Real one-loop amplitudes almost never have a trivial numerator: after the Dirac
trace and color algebra you are left with a polynomial in the loop momentum
contracted against itself and against the external momenta. This document
explains how `oneloop` represents such numerators, how it reduces them, what
is currently supported, and how a gammaloop graph numerator is translated into
the reducer's input form.

For the underlying scalar recursions this reduction rides on, see
[the reduction algorithm](02-reduction.md); for the target masters, see
[the overview](01-overview.md); for what is validated numerically, see
[the benchmarks](06-benchmarks.md).

## The numerator currency: `dot(k, q_i)`

A `oneloop` numerator is a **polynomial in Minkowski dot products** of the loop
momentum `k` with itself and with the external momenta. It is stored as a single
Symbolica `Atom` in the family:

```rust
IntegralFamily {
    // ...
    numerator: Atom,   // polynomial in oneloop::dot(k, q_i)
}
```

The building blocks are the registered symbol `oneloop::dot`, which is
**symmetric and linear in both arguments**, applied to the loop momentum
`oneloop::k` and the external momenta `oneloop::q1`, `oneloop::q2`,
`oneloop::q3`:

- `dot(k, k)` — the loop momentum squared (written `dot_ll()` internally);
- `dot(k, q_i)` — loop momentum against the `i`-th external momentum
  (written `dot_lq(i−1)` internally; 0-indexed, so `dot_lq(0)` is `dot(k, q1)`).

A rank-`r` tensor numerator is exactly a degree-`r` monomial (or a sum of them)
in these `dot` atoms. Examples that appear verbatim in the test suite:

| numerator | rank | meaning |
|---|---|---|
| `1` | 0 | scalar integral |
| `dot(k, q1)` | 1 | `k·q₁` |
| `2·dot(k,q1) − dot(k,q2) + dot(k,k)` | 2 | mixed linear + `k²` |
| `dot(k,q1) · dot(k,q1)` | 2 | `(k·q₁)²` |
| `dot(k,k)` on a box, with `[2,1,1,1]` powers | 2 + raised power | dotted box |

The scalar case is detected simply by `numerator == Atom::num(1)`; anything
else is routed to a per-topology numerator reducer.

## How a numerator reduces: dots ↔ inverse propagators

The core identity is that every reducible loop dot can be re-expressed in terms
of the **inverse propagators** `D_i` of the topology. Writing the inverse
propagator of line `i` as

```
D_i = dot(k,k) + 2·dot(k, r_i) + (r_i·r_i − m_i²),
```

where `r_i` is the external-momentum offset of line `i`, the differences
`Δ = D_i − D_0` are linear in the loop dots. Solving this linear system expresses
each `dot(k, ·)` in the reducible span as a linear combination of the `D_i` (plus
kinematic constants). In the code these are the **RSP (reducible scalar product)
rules**. For the bubble, with line 1 unshifted and line 2 shifted by `q1`
(built in `bubble_topo_general`):

```
dot(k,k) = D1 + m_a
dot(k,w) = (D2 − D1 − m_a + m_b − p²) / 2      (w = r_b − r_a)
```

and analogously for the triangle (`rule_ll`, `rule_w1`, `rule_w2`) and box
(`rule_ll`, `rule_q1..q3`). The linear solve of the Gram system
`G·α = rhs` (with `G_ij = rsp_dirs[i]·rsp_dirs[j]`) is done by **Cramer's rule
with the exact symbolic determinant** (`gram_solve`), so the Gram determinant
cancels downstream and the coefficients stay rational in `d`.

### Rank-1

For a rank-1 numerator `dot(k, q_i)`, the RSP rule immediately rewrites it as a
linear combination of `D_i`'s. Each `D_i` cancels the corresponding propagator
(lowering its exponent), producing the parent scalar integral with a shifted
power and pinched sub-topologies — which then reduce through the scalar
recursions of [the reduction algorithm](02-reduction.md). The worked example in
the MadLoop record is textbook:

```
2·dot(k,q1) = D1 − D0   ⇒   ∫ dot(k,q1)  =  ½ B0(s02;0,0) − ½ B0(s12;0,0)
```

on the massless on-shell triangle (`r1²=0`).

### Rank-2 and the transverse average

A rank-2 numerator has two loop momenta. Each is split into a **reducible part**
(the projection onto the RSP span, handled by the linear solve above) and an
**ISP part** (the component transverse to the external momenta). The reducible
parts reduce as in rank-1. The genuinely-transverse pieces are handled by a
**Passarino–Veltman transverse average**: even products of transverse loop
components are averaged over pairings of the transverse metric in
`n_t = d − |rsp|` dimensions (`isp_project`, `pairings`, `perp_metric`). Odd
transverse products integrate to zero. This keeps the coefficients rational in
`d` (the `n_t + 2j` normalizations are the only `d`-dependence introduced).

The transverse machinery is only built when a genuine ISP is actually present
(`any_isp`); for the top-level triangle and box the external momenta span the
transverse space, so numerators are purely reducible and this step is a no-op.

### The general engine

Both paths funnel through `reduce_num(topo, numerator_monos)`, which:

1. **ISP-projects** the numerator (`isp_project`);
2. substitutes the **RSP rules** (`rule_ll`, `rule_lq`) to get a polynomial in
   the `den` symbols (`oneloop::den1`, …);
3. extracts `den`-monomials (`extract_monomials`);
4. **routes** each monomial: a monomial with all `den` exponents `≤ a_i` is a
   pure propagator shift handed to the scalar reducer; a monomial that would
   raise a propagator into the numerator (`b_i < 0`) **pinches** those lines and
   carries a residual `∏ D_i^{|b_i|}` onto the surviving sub-topology
   (`pinch_to_subtopo`), which rebuilds a smaller `Topo` (box → triangle →
   bubble → tadpole) and recurses.

The `Topo` struct carries everything a level needs: the reducible directions, a
base-Gram closure over the external momenta, the RSP rules, the propagator
exponents, the per-line offsets and masses, a scalar reducer, and a pinch
routine. Sub-topologies are built by shifting into a reference line's frame
(`shift_numerator`) so the same engine handles every level uniformly.

## What is supported

Grounded in the current test suite (`reduce.rs`) and the benchmark record:

- **Scalar** numerators for any propagator count `N` (`N=1..N`).
- **Rank-1** numerators (linear in one `dot(k, q_i)` or `dot(k,k)`) — tested on
  bubble, triangle, box, pentagon, hexagon, heptagon.
- **Rank-2** numerators (degree-2 in the loop dots, e.g. `(k·q₁)²` or `dot(k,k)`)
  — tested on triangle and box, and exercised up through hexagon/heptagon in the
  reduction tests.
- **Raised propagator powers** combined with numerators (e.g. dotted box with
  exponents `[2,1,1,1]`, dotted triangle `[2,2,2]`, dotted bubble `[3,1]`).
- **On-shell massless external legs** via the `reg_delta` regularization wrapper
  (see [the reduction algorithm](02-reduction.md) and [the frontier](04-frontier.md)).

The generic cross-engine benchmark sweep went to **rank-6** (e.g. `(k²)³`,
finite and UV-divergent) and validates scalar `N=3..7`; see
[the benchmarks](06-benchmarks.md) and [`../benchmarks/README.md`](../benchmarks/README.md).
The MadLoop cross-check reproduced ~20 processes (triangle, box, pentagon,
4-gluon, loop-induced) to ~14 digits where checked. Speed of the symbolic
reduction is sub-millisecond across triangle/box/pentagon from scalar through
rank-2 (e.g. triangle rank-2 ≈ 0.27 ms, box rank-2 ≈ 0.54 ms).

### Known gaps and the higher-rank plan

- **Rank ≥ 3 on the low-point topologies** is the open generalization: the
  engine is already rank-agnostic in structure (`extract_monomials` supports
  degree up to 20 per variable, the transverse average handles any even ISP
  order), and rank-6 scalar `(k²)³` already reduces. The remaining work is
  hardening and validating the higher-rank tensor cases on triangle/box against
  oracles, and the `dot ↔ inverse-propagator` linear solve is the mechanism for
  it. See [the frontier](04-frontier.md).
- **Multiple on-shell massless legs at rank ≥ 2 on a triangle**, and massless
  internal lines — characterized and (for the on-shell-leg case) handled by the
  `reg_delta` off-shell regularization; residual edge cases are the IR frontier
  documented in [the frontier](04-frontier.md).

## The gammaloop → oneloop bridge

`bridge.rs` turns a real gammaloop one-loop graph into an `IntegralFamily` the
reducer can consume. When gammaloop is asked for a full numerator
(`save dot --output-full-numerator --do-gamma-algebra --do-color-algebra`), the
emitted numerator is a polynomial in Minkowski dot products written in
**shared-index tensor form**: a dot product `a·b` appears as
`a(mink(4,idx)) b(mink(4,idx))` (Einstein summation on the repeated index), and
`a·a` appears as a square, rather than as a `dot(...)` function. The bridge
rewrites this into `oneloop`'s `dot` convention.

The momentum map is: loop `K(0,·) → oneloop::k`; externals
`P(0/1/2,·) → oneloop::q1/q2/q3`.

### `numerator_to_dot_form`

```rust
pub fn numerator_to_dot_form(num: &Atom, heads: &GammaloopHeads) -> Atom
```

Rewrites every shared-index contraction into a `dot(...)`:

- **self-contractions** `a·a` (stored by Symbolica as `mom(mink(4,i))^2`) →
  `dot(sym, sym)`;
- **distinct contractions** `a·b` (shared-index products) → `dot(sym_a, sym_b)`;
- **metric-dot form** `g(a, b)` (which `simplify_metrics` can leave behind, for
  both self and distinct pairs) → `dot(sym_a, sym_b)`.

`GammaloopHeads` carries the four tensor heads the bridge matches against
(`loop_mom`, `external_mom`, `index`, `metric`); the real glue passes
gammalooprs's `GS.loop_mom` / `GS.external_mom` / `spenso::mink` / metric.

### `external_offset_from_lmb_rep`

```rust
pub fn external_offset_from_lmb_rep(lmb_rep: &Atom, heads: &GammaloopHeads) -> Atom
```

Extracts the external-momentum offset `r` of a propagator such that the edge
momentum is `k + r`: it zeroes every loop-momentum tensor `K(l, ·)` and maps
each external `P(j, ·)` to `q_{j+1}`. For example `K − P(0)` yields `−q1`.

### `GammaloopEdge` and `family_from_gammaloop`

```rust
pub struct GammaloopEdge { pub lmb_rep: Atom, pub mass_sq: Atom }

pub fn family_from_gammaloop(
    numerator: &Atom,
    edges: &[GammaloopEdge],
    heads: &GammaloopHeads,
) -> IntegralFamily
```

Each edge carries its loop-momentum-basis representation (`lmb_rep`) and its
mass². `family_from_gammaloop`:

- computes each edge's external **offset** via `external_offset_from_lmb_rep`;
- builds the `C(n,2)` pairwise **invariants** `(r_i − r_j)²` from those offsets
  (`invariants_from_offsets`, `square_external_momentum`);
- sets every propagator exponent to `1`;
- translates the numerator to dot form via `numerator_to_dot_form`.

This is the end-to-end path: **gammaloop graph → `IntegralFamily` → reducer →
masters**. The bridge's tests exercise it directly, e.g. a scalar massless
bubble reduces to a `B0` master, and a rank-1 bubble numerator `K·P(0)`
translates to `dot(k, q1)` and reduces through to a `B0` master plus tadpoles.

### Bridge scope

The bridge maps loop `K(0,·)` → `k` and externals `P(j,·)` → `q{j+1}`, built
**dynamically up to `MAX_MOMENTUM_ID` (8)** — so pentagons, hexagons and beyond
(up to nine-point) are handled, not only four-point topologies (extended
2026-08-21; unit-tested `maps_high_externals_for_pentagon_and_beyond` and
`invariants_handle_high_externals`, and the whole app path in
`reduce_bridge.rs` inherits it since it loops over graph edges dynamically). The
full graph extraction also depends on gammalooprs exposing the
contracted numerator atom and per-loop-edge momentum/mass publicly; the pieces
already public (`Graph::from_string`, edge mass atoms, loop-edge count) cover the
dot-export path used here. See [the app integration](05-app.md) for how this
plugs into the CLI/API.

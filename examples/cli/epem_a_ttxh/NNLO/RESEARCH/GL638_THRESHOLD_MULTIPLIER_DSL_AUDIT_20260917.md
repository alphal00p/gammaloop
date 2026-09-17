# GL638 multiplier DSL audit, 17 September 2026

Read-only audit of the current [GL638 DOT](../graphs/GL638.dot) and
[`threshold_multiplier.rs`](../../../../../crates/gammalooprs/src/integrands/process/threshold_multiplier.rs).
No multiplier, evaluator, or sampling implementation changed. The phase migration
changes the graph's global numerator factor; the six multiplier expressions are
the same real kinematic weights.

## What the literal expressions compute

At the selected variant's `star`, write `E_j=sqrt(m_j²+|q_j|²)` and
`Q=P(0,cind(0))+P(1,cind(0))`. All energies below are positive on-shell energies.
For the equal top masses in these pairs, the long rationalized components are

```
dE(i,j) = sum_k[(q_i,k-q_j,k)(q_i,k+q_j,k)] / (E_i+E_j)
        = E_i-E_j.

H_lit  = dE(4,6)+dE(3,10)
G0_lit = dE(3,10)+E13
G4_lit = dE(4,6)+E13
P_lit  = 2 E3-Q
Z_lit  = E3+E10+E13-Q

WH = H_lit² / [H_lit²+(P_lit Z_lit/Q)²]
WF = (G0_lit G4_lit)² / [(G0_lit G4_lit)²+(P_lit Z_lit)²].
```

This was checked against the actual six expression strings, not just the DOT
comments. The [standard-library script](gl638_multiplier_dsl_20260917/check_literal_expressions.py)
checks the graph's spatial momentum conservation at every vertex and evaluates
all six expressions at 43 host-shell configurations, including three soft-q13
points. With 90-digit arithmetic the largest absolute discrepancy from the short
formulas is `1.6e-89`; see [results](gl638_multiplier_dsl_20260917/results.json).
This is an algebra/notation check, not an integration or variance test.

The current centre-of-mass routing has `q12=q3`, `q5=q10`,
`q4-q6=q3-q10=q13`. The cut-preserving star satisfies
`E2+E6+E10=Q`. Consequently

```
H_lit  = (E2+E4+E12-Q)        = eta(2,4,12)
G0_lit = (E2+E6+E12+E13-Q)    = eta(2,6,12,13)
G4_lit = (E2+E4+E10+E13-Q)    = eta(2,4,10,13)
P_lit  = (E3+E12-Q)           = eta(3,12)
Z_lit  = (E3+E10+E13-Q)       = eta(3,10,13).
```

Here `eta(...)` is mathematical shorthand, not valid DSL syntax without `eset`.
In particular `2E3-Q` is justified by the equal masses and COM routing; it is not
a general edge-set identity in other external frames. Off the host shell, the
first three literal deficits equal the displayed full energy sums **minus the
host residual** `E2+E6+E10-Q`. Replacing them with a separately evaluated full
sum also changes floating-point cancellation even at a numerically solved shell.

A/U use native WH and shared 1−WH; F uses native WF and shared 1−WF. Each uses
its own star, so the values at two different off-pole stars need not sum to one.

## What the splitting achieves, and its limits

The explicit metadata belong to Cut1 `(2,6,10)`, with full parent LMB
`[3,6,7,10]`. The assignment is:

| Threshold | Shared one-loop `[7]` | Native two-loop `[3,7]` |
|---|---:|---:|
| A `(7,8)` | `1-WH` | `WH` |
| U `(8,12,14)` | `1-WH` | `WH` |
| F `(8,10,13,14)` | `1-WF` | `WF` |
| V `(3,7,14)` | absent | `1` |
| P `(3,12)` | absent | `1` |
| Z `(3,10,13)` | absent | `1` |

The other cut sides retain their natural maximal one-loop subspaces. The shared
variants join compatible cycle groups across cuts, so the common center and
subspace needed by soft LU cancellation can be used. The native variants join
P/Z/V in the same two-loop problem, preserving the common projection geometry
needed for dual cancellation at the relevant same-amplitude intersections.

In a generic soft-q13 approach, `H_lit`, `G0_lit` and `G4_lit` vanish linearly
while P/Z stay separated from zero by the massive kinematics. Consequently WH
vanishes quadratically and WF quartically: the shared one-loop variants take
over. Conversely, at a regular P or Z zero, with the respective soft marker
nonzero, the native weight tends to one and its one-loop complement vanishes
quadratically. The [earlier bias audit](../../../../../docs/research/advanced_sampling/GL638_WEIGHT_BIAS_CONTROLS.md)
details which intersections require matching residues. A global equality
`WH=WF` away from those intersections is unnecessary.

These are complementary functions at a common argument, not a claim that their
values at distinct variant stars sum to one. At the subtracted pole the
projections coincide and recover complementary leading residues. The weights
are evaluated at the radial pole and do not replace the symmetric radial
localization function used for the PV subtraction.

The construction does not remove the known hard H/Z corner. At their common
zero WH has a direction-dependent limit; its derivatives can behave as `1/R`.
Assigning the exact common point a finite value does not repair that limit.
Sampling must address the residual behavior separately. These weights alone
also do not prove bounded weights at every tangency, coincident A/P root, or
counterterm-image corner.

The phase correction rotates the full complex result. Geometry and complex-norm
cancellation checks remain applicable, but older component-specific GL638
variance studies trained the old real component, which becomes the corrected
imaginary component. They are not evidence of optimal sampling for the newly
corrected physical real part.

## Existing language and evaluation contract

| Facility | Current behavior |
|---|---|
| `eta(eset(...))`, `eta(effective|star,eset(...))` | Implemented. Default is effective; edge order is canonicalized. The edge set must select exactly one equation in the cut-owned catalogue. Missing or ambiguous external shifts are rejected. |
| `Q3(edge,index)`, `Q3(effective|star,edge,index)` | Implemented for graph-edge spatial momenta in the generation routing. Temporal component is exactly zero. Lorentz-index contractions use the Minkowski sign; a contracted spatial square is negative. |
| `P(external_position,index)` | Implemented; the first argument is the external-array position, not a graph-edge ID. |
| Model parameters, graph `params` | Implemented scalar inputs. Additional graph parameters are supplied runtime values, not local kinematic expression definitions. |
| Arithmetic, powers, `sqrt`, `if(condition,then,else)` | Inherited Symbolica expression facilities. `if` chooses `then` for nonzero condition, `else` for zero, and skips the unchosen branch. There is no separate implemented `ifzero` spelling. |
| Named kinematic definitions / aliases | No multiplier-level definitions field or user function map. The schema accepts `expression`, `symmetrize`, `opaque_derivatives`; construction supplies an empty `FunctionMap`. Existing tests called “aliases” concern equivalent Q3/eta spellings, not `let H=...`. |

The catalogue is assembled in `integrands/process/cross_section/mod.rs` from the
active cut group, related surfaces, selected threshold variants and their raised
associations. It is not the union of every foreign physical cut. Therefore the
foreign-cut equations H/G0/G4 cannot simply be substituted as eta calls. In
addition, the existing eta binder computes positive energy sums plus the stored
signed external shift directly; it does not preserve the explicit rationalized
soft differences used by these weights.
P and Z are already local thresholds, so their repeated expressions can be
written with `eta(star,eset(3,12))` and `eta(star,eset(3,10,13))`; this does not
solve the missing foreign equations or the stable-difference issue.

There is no independent CFF orientation input in this multiplier layout. Q3
uses the generation LMB's signed edge signatures; eta's signs come from its
catalogued external shift, not from applying an arbitrary orientation sign to
each positive energy. `star` is the globally reconstructed algebraic component,
including participating roots and the spectator coordinates. It does not select
another physical cut or expose a separate user-chosen CT root.

The evaluator is forced EagerOnly and rejects nonfinite or nonreal output.
The pinned Symbolica evaluator implements lazy branching through
`Instr::IfElse` and jumps over the inactive branch. Thus the actual WH expressions
assign native 0/shared 1 at their common zero without evaluating the quotient.
This defines the value at the point; it does not prove continuity there.

WF has no explicit common-zero guard. For finite physical host kinematics its
denominator is nevertheless positive: if `g=|q13|>0`, strict massive-energy
Lipschitz bounds give `G0>0` and `G4>0`; if `g=0`, `P=Z=2E10-Q` and the
three-body host bound gives `P <= -MH(2MT+MH)/Q < 0` (`-98.125 GeV` at Q=600).
This argument does not cover arbitrary off-host inputs or changed masses.

## What would shorten the expressions faithfully

These are design suggestions, not implemented syntax:

- Positive on-shell energy `E(point,edge)`, bound from the graph mass and routed
  momentum. An internal OSE atom already exists elsewhere in GammaLoop, but the
  multiplier input layout does not expose it. It is not the temporal component
  of Q3, which is zero.
- A stable energy difference `dE(point,i,j)` that retains the rationalized
  numerator, including the mass-squared difference for unequal masses. Merely
  spelling `E(i)-E(j)` shorter would lose the reason for the present formula.
- Explicit signed energy-deficit geometry at a specified effective/star point,
  with external shifts and host context, including physical-cut equations that
  are absent from the local CT catalogue. Reuse graph geometry; do not infer
  signs or silently solve a foreign cut from an edge list alone.
- Named scalar definitions or aliases for Q, H, P, Z, G0 and G4. The present
  external-energy sum is already expressible; a Q-total shorthand would need to
  define whether it means COM energy or incoming energy in the current frame.
- A quadratic partition operation only if it preserves branch semantics,
  component scaling and the supplied common-zero assignment. WH already has a
  lazy guard, so a generic safe ratio is convenience rather than missing basic
  capability. Adding a regulator or choosing a different value would change the
  prescription.

With energies, stable differences and named definitions, the literal formulas
above are already a faithful compact representation. Generalized eta support
alone would reduce some text but would not encode all these numerical contracts.

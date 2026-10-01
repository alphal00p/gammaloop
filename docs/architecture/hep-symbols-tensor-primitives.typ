= HEP symbols and tensor primitives

== Scope and result

Audit and implementation dated 2026-09-29, in the isolated
`hep-tensor-symbols` workspace, based on `ppuwtvyu/0` (`10a6c180`).
The separate `nmkklzpn` tensor-reduction migration and running notebook host
are outside this change.

All 22 Spenso-owned fields have been removed from `hep.Symbols`. Tensor
construction, patterns, factor projectors, and the canonical colour constant
belong to `symbolica.community.tensor`. A raw head remains available from
`TensorName` when an expression-level boundary is intentional.

Python uses `dirac_gamma`, `color_f`, and `color_t` for constructors, patterns,
and names. The registered symbols remain `spenso::gamma`, `spenso::f`, and
`spenso::t`.

== Complete inventory

Here `TE` means `TensorExpression`, `R` means `Representation`, and `TP`
means `TensorPattern`. Patterns remain expressions; they accept wildcards
without inferring a concrete tensor rank.

#table(
  columns: (1.1fr, 2.7fr), align: left, inset: 5pt,
  table.header([Removed `hep.Symbols` field], [Replacement]),
  [`metric`], [`TE.g(rep)(i, j)`; `TP.g(rep_pattern, i_, j_)`.],
  [`dirac_gamma`], [`TE.dirac_gamma(D)(i, j, mu)`; matching `TP.dirac_gamma`.],
  [`color_f`], [`TE.color_f(dA)(a, b, c)`; matching `TP.color_f`. Argument: adjoint dimension.],
  [`lorentz`], [`R.mink(D)`, `rep(index)`; `PortPattern.exact(rep, i_)`.],
  [`spinor`], [`R.bis(d)` and the same slot/pattern interface.],
  [`color_fundamental`], [`R.cof(Nc)`; `.dual()` selects the antifundamental representation.],
  [`color_adjoint`], [`R.coad(dA)` and its port patterns.],
  [`color_number`], [`Nc: Expression`, exported from the tensor package.],
  [`conjugate`], [`BroadcastFunction.conj()(expr)`; also accepts wildcard expressions.],
  [`dot`], [`dot(p, q)`; `TP.dot(left_, right_)`.],
  [`casimir`], [`rep.casimir(degree)`; `TP.casimir(degree_, rep_)`.],
  [`color_index`], [`rep.dynkin_index(degree)`; `TP.dynkin_index(degree_, rep_)`.],
  [`chain_in`], [`AUTO` during construction; `PortPattern.chain_in()` in contextual patterns.],
  [`chain_out`], [`AUTO` during construction; `PortPattern.chain_out()` in contextual patterns.],
  [`gamma_chain`], [`chain(left_slot, right_slot, *factors)`; `TP.chain(left_, right_, factors___)`.],
  [`levi_civita`], [`TE.levi_civita(rep, rank=4)`; `TP.levi_civita(rep_pattern, *indices)`.],
  [`charge_conjugation`], [`TE.charge_conjugation(d)`; matching `TP.charge_conjugation`.],
  [`cyclic`], [`FactorProjector.cyclic(*factors)`; `TP.cyclic(*patterns)`.],
  [`symmetric`], [`FactorProjector.symmetric(*factors)`; `TP.symmetric(*patterns)`.],
  [`antisymmetric`], [`FactorProjector.antisymmetric(*factors)`; `TP.antisymmetric(*patterns)`.],
  [`gamma_zero`], [`TE.gamma0(d)`; matching `TP.gamma0`.],
  [`trace`], [`trace(rep, *factors)` or `expr.trace()`; `TP.trace(rep_, factors___)`.],
)

`color_t` also has matching expression, pattern, and name accessors; it was
not one of the HEP fields.

== Open tensors and antisymmetry

The new expression factories `gamma0`, `charge_conjugation`, and
`levi_civita` use the existing signature/open-port constructor. It assigns
distinct unresolved port identities and retains their logical order.
Matching `TensorName` accessors return the canonical registered heads.

```python
from symbolica.community.tensor import Representation, TensorExpression, TensorName
spin = Representation.bis(4)
mink = Representation.mink(4)
C = TensorExpression.charge_conjugation(4)
epsilon = TensorExpression.levi_civita(mink)
indexed_C = C("a", "b")
indexed_epsilon = epsilon("mu", "nu", "rho", "sigma")
```

Identical arguments to an antisymmetric tensor vanish. In particular,
`TensorName.charge_conjugation()(spin, spin)` and
`TensorName.levi_civita()(mink, mink, mink, mink)` correctly produce zero.
That behavior is normal and has not been changed. Use an expression factory
for distinct unresolved ports, or pass distinct explicit slots to its name.
Repeating an explicit index also gives zero.

Epsilon rank is independent of representation dimension; arbitrary rank and
symbolic dimensions are supported. No numerical epsilon table is implied.
The existing HEP gamma-zero and charge-conjugation component data are 4D.
Symbolic factories do not create component data for other dimensions.

== Compact patterns

```python
from symbolica import S
from symbolica.community.tensor import PortPattern, TensorPattern
D_, mu_, left_, right_, factors___ = S("D_", "mu_", "left_", "right_", "factors___")
incoming = PortPattern.chain_in()
outgoing = PortPattern.chain_out()
reverse = TensorPattern.dirac_gamma(D_, outgoing, incoming, mu_)
word_pattern = TensorPattern.chain(left_, right_, reverse, factors___)
trace_pattern = TensorPattern.trace(S("rep_"), factors___)
```

Contextual ports distinguish forward and reversed channels. `AUTO` is a
concrete-construction placeholder, not an orientation for a pattern. The
trace builder retains the canonical cyclic wrapper. Casimir/Dynkin patterns
accept arbitrary representation expressions, while representation methods
remain the ordinary scalar construction API. For an unrestricted epsilon
port sequence, use `TensorPattern(TensorName.levi_civita(), ports=[factors___])`.
Compact rank-one vector expressions are accepted in contextual patterns.

== Factor projectors

```python
from symbolica.community.tensor import AUTO, FactorProjector, chain, trace
gamma = TensorExpression.dirac_gamma(4)
factors = [gamma(AUTO, AUTO, label) for label in ("mu", "nu", "rho", "sigma")]
group = FactorProjector.antisymmetric(*factors[1:3])
word = chain(spin("a"), spin("b"), factors[0], group, factors[3])
expanded = word.expand_projectors()
closed = trace(spin, FactorProjector.symmetric(*factors))
```

Groups project exactly their supplied factors; nested projector groups are
allowed. Pass a precomposed chain's individual factors explicitly. Passing
a multi-factor chain as one factor is rejected, avoiding an unintended change
to which factors are projected.

Symmetric/antisymmetric projectors average permutations with weight `1/n!`;
cyclic projectors average rotations with weight `1/n`. The Rust
`FactorProjector` owner shares coefficients with the existing
`ProjectorExpander`. `TensorExpression.expand_projectors()` retains the
interface and leaves unrelated products and sums factored. Ordinary
`expand()` continues to leave projector heads intact. Scalar multiples of a
root chain can compose and trace, including the minus sign from reversing an
antisymmetric group. Broadcast wrappers keep their separate scope.

Symbolic groups produce `TensorExpression` through `chain`/`trace`.
Component `Tensor` or `TensorNetwork` factors retain their data and produce
`TensorNetwork`. Permutations preserve spectator indices and align their
axes before addition. Stubs carry this distinction through
`FactorProjector[TensorExpression]` and `FactorProjector[TensorNetwork]`;
arbitrary mixed variadic sequences use a conservative union return type.

== Scalars, introspection, and the remaining HEP surface

`Nc` is `CS.nc` itself: a real scalar with default numerical value 3.
An unrelated `S("Nc")` lacks that identity and metadata. Generic
user-supplied colour dimensions remain supported.

Use `rep.name.name` to inspect representation names. Intentional raw-head
boundaries use `rep.to_expression().get_head()`, `rep.casimir().get_head()`,
or `TensorName.dirac_gamma().to_expression()`. These are useful for rewrite
oracles and syntax inspection.

`hep.Symbols` retains ten graph/model entries: `edge_momentum`, `denominator`,
`dimension`, `half_edge`, `polarization`, `polarization_conjugate`,
`ufo_metric`, `ufo_index`, `ufo_momentum`, and `model_conjugate`.
External and loop momenta belong to `Kinematics`; model parameters and
couplings are retrieved from their model.

== Migration and verification

The display atlas uses typed dots, chains, and traces; remaining raw examples
intentionally illustrate contextual notation. Benchmarks use the epsilon
factory/patterns, and diagnostics use representation metadata and metric
patterns. Colour-conjugation tests use checked factor groups. Charge and
gamma-zero tests use the new factories. Independent raw formulas and parser
fixtures remain expression-level oracles.

`examples/reproducers/hep_symbols_primitives.py` compares the original 24
construction, rank, and raw-head identities with canonical expression oracles.
`examples/reproducers/hep_symbols_removal_audit.py` exercises the new routes.
`crates/spynso3/tests/installed_tensor_vocabulary.py` covers antisymmetry,
scalar metadata, reversed ports, local/nested projectors, independent matrix
arithmetic, and spectator-axis alignment. The API tour includes the new
group type. Generated stubs are checked against runtime signatures and
Python typing assertions.

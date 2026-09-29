= HEP symbols and tensor primitives

== Scope and result

Audit dated 2026-09-29, based on `ppuwtvyu/0` (`10a6c180`), the parent of
`default@` when this parallel workspace was created. The separate
`nmkklzpn` reduction-API migration is outside this snapshot. Rust crate names
and symbolic `spenso::` identities remain unchanged by the Python namespace
rename to `symbolica.community.tensor`. The follow-up API rename in this
parallel branch uses `dirac_gamma`, `color_f`, and `color_t` for constructors, patterns,
and names; the registered symbols remain `spenso::gamma`, `spenso::f`, and `spenso::t`.

`crates/feynkit-py/src/symbols.rs` exposes 22 Spenso-owned heads through
`hep.Symbols`, each as an ordinary Symbolica `Expression`. Most concrete
construction already belongs to the tensor API. A raw head and a typed
constructor are different interfaces: a constructor knows dimensions, slot
order, duality, and the tensor result type; a raw head also permits patterns
and incomplete syntax. Replacing every accessor with a callable would hide
that distinction rather than remove it.

The audit records an executable comparison; the follow-up changes the API
spellings above without changing the running notebook environment. The comparison is
`examples/reproducers/hep_symbols_primitives.py`. It checks 24 exact-expression,
rank, and raw-head equalities in a combined community host. Its
`to_expression()` calls are the oracle boundary; the typed examples themselves
do not require conversion.

== Complete inventory

In this table, `TE` means `TensorExpression`, `R` means `Representation`,
`rep` is a representation, and `lib` is `TensorLibrary.hep_lib_atom()`.
Quoted library keys are fully qualified because a bare string uses Symbolica's
current namespace.

#table(
  columns: (1.15fr, 2.4fr, 1.6fr),
  align: left,
  inset: 5pt,
  table.header([`hep.Symbols`], [Existing tensor primitive], [Boundary or gap]),
  [`metric`], [`TE.g(rep)(i, j)`], [Typed factory; accepts an explicit dual partner.],
  [`dirac_gamma`], [`TE.dirac_gamma(D)(i, j, mu)`], [Bispinor, bispinor, Lorentz order.],
  [`color_f`], [`TE.color_f(dA)(a, b, c)`], [Argument is adjoint dimension, not number of colours.],
  [`lorentz`], [`R.mink(D)` and `rep(index)`], [Representation and slot.],
  [`spinor`], [`R.bis(d)` and `rep(index)`], [Representation and slot.],
  [`color_fundamental`], [`R.cof(Nc)` and `rep(index)`], [Use `rep.dual()` for antifundamental slots.],
  [`color_adjoint`], [`R.coad(dA)` and `rep(index)`], [Representation and slot.],
  [`color_number`], [`rep.dimension` for a chosen colour representation], [No tensor accessor for the canonical `spenso::Nc` constant.],
  [`conjugate`], [`BroadcastFunction.conj()(expr)`], [Preserves a symbolic tensor's interface.],
  [`dot`], [`dot(p, q)` or `p.dot(q)`], [Rank-one tensor operands yield a tensor scalar.],
  [`casimir`], [`rep.casimir(degree)`], [Returns a scalar `Expression`; default degree is 2.],
  [`color_index`], [`rep.dynkin_index(degree)`], [Returns a scalar `Expression`; default degree is 2.],
  [`chain_in`], [`AUTO` ports supplied to `chain(...)`], [Builder emits the internal wiring marker.],
  [`chain_out`], [`AUTO` ports supplied to `chain(...)`], [Builder emits the internal wiring marker.],
  [`gamma_chain`], [`chain(left_slot, right_slot, *factors)`], [Typed builder; retains order and spectator ports.],
  [`levi_civita`], [`TensorName("spenso::epsilon")(*distinct_slots)`], [Works when explicitly indexed; dedicated factory missing.],
  [`charge_conjugation`], [`lib["spenso::charge_conjugation"].expression()`], [Typed 4D reference; dedicated factory missing.],
  [`cyclic`], [`trace(rep, *factors)` for ordinary traces], [Builder inserts it; standalone cyclic projector lacks a typed builder.],
  [`symmetric`], [Rust `spenso::shadowing::projectors::sym`], [Typed Python factor-projector constructor missing.],
  [`antisymmetric`], [Rust `spenso::shadowing::projectors::antisym`], [Typed Python factor-projector constructor missing.],
  [`gamma_zero`], [`lib["spenso::gamma0"].expression()`], [Typed 4D reference; dedicated factory missing.],
  [`trace`], [`trace(rep, *factors)` or `expr.trace()`], [Construct a trace word or close a tensor channel.],
)

== Idiomatic construction

```python
from symbolica.community.tensor import (
    AUTO, BroadcastFunction, Representation, TensorExpression,
    TensorLibrary, TensorName, chain, dot, trace,
)

mink = Representation.mink(4)
spin = Representation.bis(4)
fund = Representation.cof(3)

gamma = TensorExpression.dirac_gamma(4)
g = TensorExpression.g(mink)
p = TensorName.vector("example::p")(mink)
q = TensorName.vector("example::q")(mink)
pq = dot(p, q)
conjugate_p = BroadcastFunction.conj()(p)
cf = fund.casimir()
tr = fund.dynkin_index()

factors = (gamma(AUTO, AUTO, "mu"), gamma(AUTO, AUTO, "nu"))
word = chain(spin("a"), spin("b"), *factors)
closed_word = trace(spin, *factors)

library = TensorLibrary.hep_lib_atom()
C = library["spenso::charge_conjugation"].expression()
gamma0 = library["spenso::gamma0"].expression()
indexed_C = C("a", "b")
```

`BroadcastFunction.conj()` is componentwise conjugation. It is not a synonym
for `dirac_adjoint()`, which can reverse matrix chains and apply the
representation-aware adjoint. Nor should callers assume that inherited
Symbolica `conj()` preserves a `TensorExpression` result.

`trace(...)` leaves any uncontracted Lorentz ports open: tracing the two gamma
factors above yields rank two. It does not evaluate the gamma trace.

== Missing constructors and unsafe shortcuts

The useful convenience additions belong on `TensorExpression`, alongside its
existing factories: `gamma0`, `charge_conjugation`, and `epsilon`. Their
implementations should use `AGS.gamma0_strct`, `AGS.charge_conjugation`, and
`EPSILON_SYMBOL`/the existing epsilon builder. They must assign distinct
unresolved port identities before Symbolica normalization, retain logical
port order, and validate the representation requirements. The existing Rust
gamma-zero structure builder already supplies its canonical head and shape.
The HEP library route above selects its actual stored four-dimensional data;
it is not a general symbolic-dimension factory.

The generic shortcut is unsafe for antisymmetric tensors with repeated open
representations:

```python
TensorName("spenso::epsilon")(mink, mink, mink, mink)
TensorName("spenso::charge_conjugation")(spin, spin)
```

Both produce an underlying zero in the inspected runtime. Their argument
lists contain identical unindexed representations when the antisymmetric
normalizer runs. Explicit, distinct slots work. A factory must avoid this
premature collapse; merely moving the raw head from HEP to another class is
insufficient. No test here asserts that this faulty behavior should persist.

The factor projectors `sym`, `antisym`, and `cyclic` already have canonical
Rust constructors in `crates/spenso/src/shadowing/projectors.rs`. They are
normalized projectors over ordered matrix factors. `TensorName(...,
is_symmetric=True)` declares an atomic tensor's symmetry; it is not a
replacement for a projector over several matrix factors. Python needs a
factor-projector API only for callers that explicitly request these objects;
ordinary chain and trace construction already handles its wiring wrappers.

`color_number` needs a separate decision about scalar conventions. A
caller-supplied `Nc` and `Representation.cof(Nc).dimension` suffice for generic
SU(N) calculations. They do not retrieve the canonical `CS.nc` symbol:
`spenso::Nc` also carries reality and default numerical evaluation metadata
(value 3). Replacing it with an arbitrary `S("Nc")` would change identity and
numerical behavior. Expose the canonical constant from the tensor/color owner
when that identity is required.

== Patterns, introspection, and ownership

For rewriting, use `TensorPattern.g`, `TensorPattern.dirac_gamma`,
`TensorPattern.color_f`, `TensorPattern.vector`, and `PortPattern`. For raw head
introspection, `TensorName.g().to_expression()`, `TensorName.dirac_gamma().to_expression()`,
and `TensorName.color_f().to_expression()` already reproduce the HEP handles.
The predefined `TensorName` accessors are intentionally not concrete callable
factories. Representation identity is available through `rep.name`; do not
construct a raw `mink(...)` call merely to inspect a port's representation.

`Symbols.edge_momentum`, `denominator`, `dimension`, `half_edge`,
`polarization`, `polarization_conjugate`, and the four `ufo_*`/model-conjugation
entries are graph or model vocabulary. They are not among the 22
Spenso-owned entries audited here. External and loop momenta already belong
to `Kinematics` in this baseline.

Recommended migration order:

+ Migrate concrete construction to the existing typed primitives and raw
  rewrites to the pattern API. Keep scalar invariants scalar.
+ Add the missing safe tensor factories in the tensor package, plus matching
  pattern shortcuts where real callers need them.
+ Resolve canonical colour constants and explicit factor projectors at their
  tensor/color owners.
+ Remove Spenso-owned exports from `hep.Symbols` and update generated stubs,
  documentation, and callers together. Do not create a second tensor-symbol
  registry in HEP or keep compatibility wrappers for construction.

= HEP symbols and tensor primitives

== Scope and result

Audit dated 2026-09-29, based on `ppuwtvyu/0` (`10a6c180`), the parent of
`default@` when this parallel workspace was created. The separate
`nmkklzpn` reduction-API migration is outside this snapshot. Rust crate names
and symbolic `spenso::` identities remain unchanged by the Python namespace
rename to `symbolica.community.tensor`. The follow-up API rename in this
parallel branch uses `dirac_gamma`, `color_f`, and `color_t` for constructors, patterns,
and names; the registered symbols remain `spenso::gamma`, `spenso::f`, and `spenso::t`.
The removal-readiness follow-up below was checked against `zpmpryrx`
(`a70421b6`) and its compiled Python host. The concurrent `default@`
(`43ba8c7efb90` at inspection) still lacks the additional factories and pattern
builders identified here. No runtime API is changed by this follow-up audit.

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

== Removal readiness: complete surface

Deleting aliases is easy; retaining a usable API for both concrete tensors and
rewrite rules needs a broader check. The original 24 comparisons establish
construction equivalence. `examples/reproducers/hep_symbols_removal_audit.py`
adds wildcard matching checks and reports the current gaps without making the
faulty behaviors permanent regression expectations.

The following table accounts for all 22 entries. “Ready” means an existing
semantic API covers construction and the ordinary pattern/introspection use;
advanced pattern gaps are stated explicitly. All proposed names below are
design recommendations, not APIs that already exist.

#table(
  columns: (1.3fr, 2.7fr),
  align: left,
  inset: 5pt,
  table.header([Entries in `hep.Symbols`], [Removal requirement]),
  [`metric`, `dirac_gamma`, `color_f` (3)],
  [Ready: matching `TensorExpression`, `TensorPattern`, and `TensorName` factories. Generic `TensorPattern(name, ports=[ports___])` covers arbitrary port sequences.],
  [`lorentz`, `spinor`, `color_fundamental`, `color_adjoint` (4)],
  [Ready: `Representation.mink/bis/cof/coad`, `PortPattern.exact`, and tagged `PortPattern.any/dualizable`. Wildcard dimensions already work.],
  [`conjugate` (1)],
  [Ready: `BroadcastFunction.conj()` works on tensors and ordinary wildcard expressions. Its explicit head conversion preserves the registered symbol.],
  [`casimir`, `color_index` (2)],
  [Ready for a specified representation family: `rep.casimir(k_)` and `rep.dynkin_index(k_)`, including wildcard dimensions. Add pattern shortcuts for matching an arbitrary representation expression.],
  [`dot` (1)],
  [Concrete `dot(p, q)` already works. Add a compact-dot pattern constructor accepting arbitrary operand patterns.],
  [`gamma_chain`, `trace` (2)],
  [Concrete `chain` and `trace` already work. Add pattern constructors accepting wildcard endpoints/representations and factor sequence wildcards; preserve canonical cyclic trace storage.],
  [`chain_in`, `chain_out` (2)],
  [Concrete `AUTO` ports already replace direct marker use. Pattern builders also need contextual input/output ports, including transposed factors. No public HEP marker aliases are needed.],
  [`gamma_zero`, `charge_conjugation`, `levi_civita` (3)],
  [Add safe named expression factories, matching pattern factories, and `TensorName` accessors. The latter two require a construction fix for repeated unresolved ports.],
  [`symmetric`, `antisymmetric`, `cyclic` (3)],
  [Add checked factor-projector construction and matching patterns. Ordinary traces already supply their cyclic wrapper, but explicit symmetric/antisymmetric words need their own public route.],
  [`color_number` (1)],
  [Export the existing canonical scalar `Nc` from the tensor package, retaining its identity, reality attribute, and default numerical evaluation.],
)

== Recommended additions

=== Three tensor families

Add corresponding factories to `TensorExpression`, `TensorPattern`, and
`TensorName` for `gamma0`, `charge_conjugation`, and `levi_civita`. The proposed
Python `levi_civita` spelling still denotes `spenso::epsilon`.

Suggested construction signatures are `TensorExpression.gamma0(4)`,
`TensorExpression.charge_conjugation(4)`, and
`TensorExpression.levi_civita(rep, rank=4)`. Keep epsilon rank separate from
the representation dimension: the existing Rust epsilon builder supports
arbitrary rank, including symbolic representation dimensions. Its numerical
component convention must remain explicit; this audit does not propose
populating a new epsilon data table or changing the existing HEP convention.
The stored gamma-zero and charge-conjugation component matrices are 4D; a
symbolic factory must not imply that numerical data or simplification
identities exist for every dimension.

Both the underlying construction and the named factories must use the shared
signature/open-port path before antisymmetric normalization. Today:

```python
spin = Representation.bis(4)
mink = Representation.mink(4)
TensorName("spenso::charge_conjugation")(spin, spin)  # rank 2, value 0
TensorName("spenso::epsilon")(mink, mink, mink, mink)  # rank 4, value 0
```

The loss occurs while forming the atom with identical unresolved arguments;
assigning distinct indices afterwards cannot recover it. Distinct explicit
slots work, and the library's charge-conjugation signature works. Fix the
shared generic path as well as exposing the safe named factories. Tests should
compare explicit indexing and reindexing, antisymmetry, port order, and the
existing library components without blessing the current zero result.

=== Compact and scalar pattern builders

Extend `TensorPattern` with `dot`, `chain`, and `trace` builders. They must
accept ordinary and sequence wildcards without asking the concrete tensor
inference engine to assign them a rank. Current concrete calls
`dot(x_, y_)`, `chain(start, end, factors___)`, and
`trace(rep, factors___)` all reject those unrestricted patterns.

The chain/trace pattern interface must cover selected factors as well as a
whole factor sequence, reversed matrix channels, and compact vector insertions.
It must lower contextual input/output ports to the canonical markers internally.
Simply teaching patterns to use `AUTO` is insufficient unless the builder can
also distinguish the reversed orientation; currently even
`TensorPattern.dirac_gamma(4, AUTO, AUTO, mu_)` is rejected.

Add pattern-only `casimir(degree_pattern, representation_pattern)` and
`dynkin_index(degree_pattern, representation_pattern)` alongside these
builders. They close the unrestricted `cas(k_, rep_)` / `idx(k_, rep_)` case.
No new scalar construction API is necessary: representation methods already
accept wildcard degrees and dimensions, and the result should remain an
ordinary scalar expression.

No additional conjugation builder is necessary. For example,
`BroadcastFunction.conj()(x_)` already constructs the desired pattern.

=== Factor projectors

Expose the existing normalized Rust factor projectors through a checked
`FactorProjector` value with `symmetric`, `antisymmetric`, and `cyclic`
constructors. Let existing `chain` and `trace` accept it as a factor group.
This preserves partial projections such as `chain(left, right, A,
FactorProjector.antisymmetric(B, C), D)`; a single symmetry flag on the whole
chain would not cover that case.

A projector group carries ordered factors and a compatible matrix channel;
its contextual `in`/`out` syntax is not a standalone tensor. In the current
host, wrapping `sym(gamma(in,out,mu), gamma(in,out,nu))` in `TensorExpression`
fails, while placing the same group inside a complete trace succeeds. The
binding should construct and validate the complete chain/trace through the
existing Rust owners rather than require a raw-expression detour.

Add matching pattern constructors, and expose the existing Rust
`ProjectorExpander` as a tensor-preserving `expand_projectors()` operation.
Ordinary `expand()` currently leaves the compact projector intact. Explicit
expansion must retain the existing normalized coefficients: `1/n!` for the
full symmetric/antisymmetric sums and `1/n` for cyclic averaging. Keep unrelated
products and sums factored. Pattern values remain expressions; complete
symbolic words remain `TensorExpression`, and component/network operands
should follow the existing `TensorNetwork` dispatch policy.

=== Canonical colour constant

Export `Nc: Expression` directly from `symbolica.community.tensor`, using the
existing `CS.nc`. It is a scalar constant, so neither a tensor factory nor a
second `Symbols` class is appropriate. In the inspected runtime,
`Symbols.color_number.evaluate({})` returns `3+0j`, and the head is real.
`S("Nc")` is a different namespaced symbol and lacks that evaluation metadata.
Keep user-supplied generic colour dimensions supported independently.

== Existing APIs to reuse

There is no need to add four representation-symbol accessors. Port patterns
can already be made with `PortPattern.exact(Representation.mink(D_), i_)`.
For name comparisons, use `rep.name.name`; for an intentional raw-head
boundary, `rep.to_expression().get_head()` already returns the exact registered
symbol. Likewise, `rep.casimir().get_head()` returns the Casimir head.
`RepresentationName.to_expression()` does not exist, but it is not a removal
blocker given these existing routes.

Typed multiplication and contraction already construct compact slash factors
for the display atlas. In this snapshot, the simplification result is explicitly
materialized to a `TensorExpression` as follows:

```python
slash = (
    TensorExpression.dirac_gamma(4)(AUTO, AUTO, "mu") * p("mu")
).contract().to_expression()
word = chain(spin("i"), spin("j"), slash)
```

This is equivalent to the compact gamma insertion built with raw `chain`,
`in`, and `out` heads. Direct `gamma(AUTO, AUTO, p)` is currently rejected,
but adding that shortcut is not necessary to remove `Symbols`.

== Migration sites and acceptance

- `examples/notebooks/spenso_idenso_display.py`: replace raw dot/chain/trace
  construction and endpoint markers with the existing builders; use the new
  gamma-zero factory. Keep raw syntax only where it is intentionally the object
  being demonstrated.
- `examples/notebooks/tensor_benchmark_cases.py`: use the safe epsilon signature
  for component registration and its pattern factory for leaf replacement.
- `examples/notebooks/tensor_benchmark_diagnostics.py`: use representation
  metadata and metric patterns/head introspection instead of HEP aliases.
- `examples/notebooks/feynkit/07_tensor_reduction_marimo.py`: replace the raw
  Lorentz head with representation/port constructors during its API migration.
- `crates/feynkit-py/tests/installed_color_conjugation.py`: replace raw
  cyclic-to-symmetric/antisymmetric rewrites with the factor-projector API.
  `installed_charge_conjugation.py` and the FeynCalc fixtures also contain direct
  `S("spenso::...")` calls and must be audited, not merely grepped for `Symbols`.

Before deletion, verify all 22 replacements in an installed tensor-only host,
including wildcard matching and head metadata; cover reversed chains, cyclic
traces, partial projections, and open antisymmetric tensors. The combined host
used here registers both packages, so its successful comparisons alone do not
certify import-order independence. Regenerate stubs/docs and remove the HEP
fields and caller imports in the same implementation. The final `hep.Symbols`
surface retains the ten graph/model entries listed above. `color_t` needs no
HEP alias: its three tensor factories already exist, and it was never one of
the 22 HEP fields.

# Spenso Python API

Spenso is available from `symbolica.community.spenso`. Its central Python types
separate concrete tensor structure, pattern syntax, stored data, and execution:

- `TensorExpression` carries a Symbolica expression and its ordered external
  tensor interface.
- `TensorPattern` and `PortPattern` build tagged Symbolica patterns without
  requiring a concrete tensor interface.
- `Tensor` owns dense or sparse data and has a `TensorExpression` as its exact
  structure.
- `TensorNetwork` represents operations in which concrete tensors or existing
  networks participate.

## Symbolic expressions

```python
from symbolica.community.spenso import _, Representation, TensorExpression, TensorName

mink = Representation.mink(4)
bis = Representation.bis(4)

p = TensorName.vector("p", is_linear=True, tags=["kinematics"])
gamma = TensorExpression.gamma(4)

square = p(1, mink) * p(1, mink)
line = gamma("mu", _, _) * gamma("nu", _, _)
dirac_trace = line.trace()

print(square.format_tensor())
print(dirac_trace.to_typst())
display(dirac_trace.formatted())
```

Calling a user-defined `TensorName` returns a `TensorExpression`. Scalar key
arguments come first and structural `Slot` or `Representation` arguments come
after them. A `Slot` is an explicit port; a `Representation` is an unresolved
port. They may be mixed in one call:

```python
mu = mink("mu")
mixed = TensorName("A")(7, mu, bis)

assert mixed.rank == 2
assert mixed.interface == (mu, bis)

# TensorExpression indexing addresses unresolved ports only.
indexed = mixed(bis("i"))
```

`_` (also exported as `AUTO`) leaves that local unresolved port open. It is
construction syntax, not a shared Einstein index:

```python
partly_indexed = gamma(mu, _, _)
fully_indexed = partly_indexed(bis("i"), bis("j"))
```

`TensorName.vector(...)` enforces one final structural argument. General tensor
names may have any number of structural ports, including none:

```python
mass = TensorName("mass")(1)       # rank-zero TensorExpression
momentum = p(1, mink)              # rank-one TensorExpression
```

### Predefined tensors

Predefined tensors have typed factories on `TensorExpression`. A factory fixes
the representations and returns an expression whose ports are unresolved;
calling that result supplies indices in logical interface order:

```python
generator = TensorExpression.t(8, 3)
structure_constant = TensorExpression.f(8)

T_aij = generator("a", "i", "j")
f_abc = structure_constant("a", "b", "c")
gamma_muij = TensorExpression.gamma(4)("mu", "i", "j")
```

The available factories and their logical interfaces are:

| Factory | Logical interface |
| --- | --- |
| `TensorExpression.g(rep)` | `rep, rep` |
| `TensorExpression.flat(rep)` | `rep, rep` |
| `TensorExpression.gamma(D)` | `mink(D), bis(4), bis(4)` |
| `TensorExpression.gamma5(D)` | `bis(D), bis(D)` |
| `TensorExpression.projm(D)` | `bis(D), bis(D)` |
| `TensorExpression.projp(D)` | `bis(D), bis(D)` |
| `TensorExpression.sigma(D)` | `mink(D), mink(D), bis(4), bis(4)` |
| `TensorExpression.f(DA)` | `coad(DA), coad(DA), coad(DA)` |
| `TensorExpression.t(DA, DF)` | `coad(DA), cof(DF), coaf(DF)` |

Factory dimensions belong to these representations. They are not scalar tensor
arguments and therefore do not add fields to tensor-library keys. The
corresponding `TensorName.g()`, `TensorName.gamma()`, and other predefined name
accessors expose the raw fixed heads for inspection and pattern construction;
calling those names directly as concrete tensors is an error.

## Tensor patterns

`TensorPattern` and `PortPattern` construct matching syntax rather than concrete
tensors. They inherit Symbolica expression behavior, so they can be passed
directly to `Expression.replace`, tensor-network replacement, and conditions.
Pattern dimensions, indices, and variadic tails remain ordinary Symbolica
wildcards:

```python
from symbolica import S
from symbolica.community.spenso import PortPattern, TensorName, TensorPattern

D_, mu_, i_, j_, args___ = S("D_", "mu_", "i_", "j_", "args___")

gamma_pattern = TensorPattern.gamma(D_, mu_, i_, j_)
fixed_head = TensorPattern(
    TensorName("A"),
    args=[args___],
    ports=[PortPattern.self_dual("R_", D_, i_)],
)
any_tensor = TensorPattern.any("T_", ports=[args___])
any_vector = TensorPattern.vector(
    "P_", ports=[PortPattern.any("R_", D_, i_)]
)
```

The general constructor always emits scalar `args` before structural `ports`.
A port may be a concrete `Representation` or `Slot`, a `PortPattern`, or an
arbitrary expression such as a sequence wildcard. `TensorPattern.any` constrains
its head to tensor-tagged symbols, while `TensorPattern.vector` additionally
requires the rank-one tag. Wildcard head names end in exactly one underscore.
Dimensions, indices, and scalar arguments use numbers or Symbolica expressions;
bare strings are reserved for the wildcard-head names rather than implicitly
creating pattern variables.

Use `PortPattern.exact(rep, index=None)` for a fixed representation.
`PortPattern.any`, `.self_dual`, and `.dualizable` constrain a wildcard
representation head by its Spenso tags. Passing `dual=True` to `.dualizable`
matches the dual orientation. Omitting `index` in any of these methods builds a
stripped representation pattern.

Every predefined tensor also has a pattern shortcut whose indices follow the
same logical ordering as its concrete factory:

```python
TensorPattern.g(rep_pattern, i_, j_)
TensorPattern.flat(rep_pattern, i_, j_)
TensorPattern.gamma(D_, mu_, i_, j_)
TensorPattern.gamma5(D_, i_, j_)
TensorPattern.projm(D_, i_, j_)
TensorPattern.projp(D_, i_, j_)
TensorPattern.sigma(D_, mu_, nu_, i_, j_)
TensorPattern.f(DA_, a_, b_, c_)
TensorPattern.t(DA_, DF_, a_, i_, j_)
```

The builders translate logical ordering to Spenso's canonical atom ordering.
Callers therefore do not need to encode internal storage permutations in a
replacement rule.

## Tensor-aware algebra

Addition and subtraction align fully indexed interfaces by abstract index, so
permuted terms such as `A(i, j) + A(j, i)` are valid. Interfaces containing
unresolved ports must have the same logical order. Scalar zero is an additive
identity, including the initial zero used by Python's `sum()`; other scalar
addends remain invalid for non-scalar tensors. Multiplication uses the interface
to choose one unambiguous operation:

1. Compose unique two-ended propagation channels into an ordered chain.
2. Otherwise contract the unique compatible port pair (rank one by rank one is
   a dot product).
3. Form an outer product when no compatible pair exists.
4. Report every candidate position when the choice is ambiguous.

The explicit forms use positions in `TensorExpression.interface`:

```python
from symbolica.community.spenso import chain, dot, trace

a = TensorName("A")(mink, bis, bis)
b = TensorName("B")(mink, bis, bis)

outer = a.outer(b)
one_link = a.contract(b, left=0, right=0)
matrix_product = a.compose(b, left=(1, 2), right=(1, 2))
closed = matrix_product.trace(channel=(0, 1))

scalar = dot(p(1, mink), p(2, mink))
open_line = chain(bis("i"), bis("j"), gamma(mu, _, _), gamma(mu, _, _))
closed_line = trace(bis, gamma(mu, _, _), gamma(mu, _, _))
```

Tensor-derived scalars remain `TensorExpression` instances with an empty
interface. This retains tensor-aware dispatch and semantic display.

`expand()` is a structure-preserving operation: it expands the Symbolica atom,
validates every summand, and returns a `TensorExpression` with the same ordered
interface and name metadata. Use `to_expression()` when an operation should
deliberately leave tensor-aware dispatch, and `as_tensor()` or `reinfer()` to
validate Spenso syntax and restore it:

```python
from symbolica import S
from symbolica.community.spenso import as_tensor

x = S("x")
expanded = ((1 + x) * p(1, mink)).expand()

plain = expanded.to_expression()
restored = as_tensor(plain)
```

## Dense and sparse data

`Tensor.dense` and `Tensor.sparse` accept a named `TensorExpression`. Its
interface may be unresolved, explicit, or mixed. Input data is row-major in the
logical order shown by `expression.interface`; Spenso's internal canonical axis
permutations do not change that public layout.

```python
from symbolica.community.spenso import Tensor

euc = Representation.euc(2)
matrix_name = TensorName("matrix")

open_matrix = matrix_name(euc, euc)
dense = Tensor.dense(open_matrix, [1.0, 2.0, 3.0, 4.0])

mixed_matrix = matrix_name(euc("row"), euc)
sparse = Tensor.sparse(mixed_matrix, float)
sparse[[0, 1]] = 3.0

assert dense.structure().interface == open_matrix.interface
print(dense[[1, 0]])               # concrete data access
```

`TensorExpression.__getitem__` maps flat logical positions to coordinates (and
coordinates back to a flat position). `Tensor.__getitem__` accesses concrete
data. Abstract indexing uses `()` or `.index()`:

```python
network = dense(euc("mu"), euc("nu"))
```

Naming a composite expression attaches data identity metadata without rewriting
its Symbolica atom. A parameterized rank-zero tensor expression can supply both
the name and scalar key arguments:

```python
left = TensorName.vector("left")(euc)
right = TensorName.vector("right")(euc)
K = TensorName("K")

definition = left.outer(right).with_name(K(1))
composite = Tensor.dense(definition, [0.0, 1.0, 1.0, 0.0])
```

Public constructors require a name. A computed network result may remain
unnamed; call `tensor.with_name(...)` only when it will later be registered.

## Automatic networks

Operations involving only tensor expressions return `TensorExpression`. As soon
as a `Tensor` or `TensorNetwork` participates, the same operation and ambiguity
rules return `TensorNetwork`. Execution remains explicit:

```python
u = Tensor.dense(TensorName.vector("u")(euc), [1.0, 2.0])
v = Tensor.dense(TensorName.vector("v")(euc), [3.0, 4.0])

network = u * v                  # inferred dot-product network
network.execute()
result = network.result_tensor() # unnamed rank-zero Tensor is valid
```

The same promotion applies to `outer`, `contract`, `compose`, `trace`, `dot`,
`chain`, and broadcast functions. A symbolic expression can be parsed directly:

```python
network = dirac_trace.to_network()
network.execute()
```

An optional `TensorLibrary` resolves symbolic tensor leaves while parsing.

## Tensor libraries

Libraries register ordinary `Tensor` objects. Registration requires a named,
fully unresolved interface (or a scalar). A key contains the tensor symbol,
its scalar arguments, and its canonical representations, so parameterized
families remain distinct:

```python
from symbolica.community.spenso import TensorLibrary

momenta = TensorLibrary()
momenta.register(Tensor.dense(p(1, mink), [1.0, 0.0, 0.0, 1.0]))
momenta.register(Tensor.dense(p(2, mink), [0.0, 1.0, 1.0, 0.0]))

p1 = momenta[p(1, mink)]         # precise lookup
p2 = momenta[p(2, mink)]
```

`library[p]` is a convenience only when exactly one registered key has that
symbol; otherwise the error lists the available signatures.

## Elementwise functions

`BroadcastFunction` is separate from tensor symbols. Symbolic broadcasting
needs no registration. A concrete network uses a callback in
`TensorFunctionLibrary`:

```python
from symbolica.community.spenso import BroadcastFunction, TensorFunctionLibrary

sqrt = BroadcastFunction("sqrt", is_real=True, tags=["kinematics"])
functions = TensorFunctionLibrary()
functions.register(sqrt, lambda value: value**0.5)

values = Tensor.dense(TensorName.vector("values")(euc), [1.0, 4.0])
network = sqrt(values)
network.execute(function_library=functions)
root_values = network.result_tensor()
```

`BroadcastFunction.conj()` returns the registered complex-conjugation builtin.
An unknown concrete callback produces a targeted execution error.

## Display

`TensorExpression`, `Tensor`, and `TensorNetwork` share semantic display:

```python
text = dirac_trace.format_tensor(show_dimensions=False)
typst = dirac_trace.to_typst(show_dimensions=False)
rich = dirac_trace.formatted(show_dimensions=False)
```

Module-level `format_tensor`, `to_typst`, and `formatted` functions also accept
ordinary Symbolica expressions containing valid Spenso syntax. Dimensions are
shown only when requested. Dots are infix, chain order is preserved, traces use
`Tr`, and internal chain wiring labels are not exposed.

## Symbolic simplification

Community simplifiers operate on ordinary Symbolica expressions. Cross the
boundary explicitly, then reinfer when the result is still valid tensor syntax:

```python
from symbolica.community.idenso import simplify_metrics

simplified = simplify_metrics(dirac_trace.to_expression())
tensor_simplified = as_tensor(simplified)
```

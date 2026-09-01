# Spenso Python API

Spenso is available from `symbolica.community.spenso`. Its public Python API has
one symbolic structured type, one data-bearing tensor type, and one executable
network type:

- `TensorExpression` carries a Symbolica expression and its ordered external
  tensor interface.
- `Tensor` owns dense or sparse data and has a `TensorExpression` as its exact
  structure.
- `TensorNetwork` represents operations in which concrete tensors or existing
  networks participate.

## Packaging and kernel ownership

`spynso3` is a component of the combined
[Symbolica Community](https://github.com/symbolica-dev/symbolica-community) Python
extension, not an independent native extension. The combined build creates
`symbolica.core` and registers `SpensoModule` as
`symbolica.community.spenso_native` in that same shared library. This is
required: Symbolica expressions, interned symbols, and their PyO3 classes must
all belong to one Symbolica kernel and one Python `Expression` type.

GammaLoop owns the Spynso source and crate-local copies of Tydenso's bundled
`render.typ` and `notation.typ` assets. Tydenso's check verifies that the two
copies remain byte-identical. Symbolica Community owns the wheel that exposes them as
`symbolica.community.spenso`. Its Cargo graph must resolve exactly one
Symbolica source revision for the core, Spynso, Spenso, and Idenso. Do not link
Spynso into `gammaloop._gammaloop` or publish it as a second native library;
either would create a separate Symbolica kernel.

The public Python package is a small wrapper around the registered native
submodule. The host adds `initialize_module` before Spynso registers its public
API, so that initializer must remain in the native export list used here:

```python
from ..spenso_native import *

initialize_module()
```

The `gammaloop[typst-display]` extra only supplies the optional Typst compiler
used by `to_html` and `to_svg`; it does not install a second Spynso binary.

## Symbolic expressions

```python
from symbolica.community.spenso import _, Representation, TensorName

mink = Representation.mink(4)
bis = Representation.bis(4)

p = TensorName.vector("p", is_linear=True, tags=["kinematics"])
gamma = TensorName.gamma()(mink, bis, bis)

square = p(1, mink) * p(1, mink)
line = gamma.index(mink("mu"), _, _) * gamma.index(mink("nu"), _, _)
dirac_trace = line.trace()

print(square.format_tensor())
print(dirac_trace.to_typst())
display(dirac_trace.formatted())
```

Calling a `TensorName` always returns a `TensorExpression`. Scalar key
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

## Tensor-aware algebra

Addition and subtraction require identical ordered interfaces. Multiplication
uses the interface to choose one unambiguous operation:

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

`DisplaySettings` selects the ports, Schoonschip, or call tensor layout and
controls dimensions, parentheses, commas, symbol scripts, and index/factor
spacing. The original positional Boolean remains supported, so
`to_typst(True)` is equivalent to requesting dimensions explicitly:

```python
from symbolica.community.spenso import DisplaySettings

source = dirac_trace.to_typst(
    settings=DisplaySettings(show_dimensions=True, parentheses=False)
)
```

`to_typst` and `format_tensor` remain source-formatting APIs and emit the
existing ports layout. They accept its native printing controls, but reject
Schoonschip, call, and custom-spacing settings with guidance to use the Typst
renderer. Those layouts depend on `notation.typ`, so use them through HTML,
SVG, or rich notebook display:

Install the optional renderer to compile the same semantic display to HTML or
SVG:

```bash
pip install 'gammaloop[typst-display]'
```

```python
compact = DisplaySettings.schoonschip()
html = dirac_trace.to_html(settings=compact)
svg = dirac_trace.to_svg(settings=compact)
rich = dirac_trace.formatted(settings=compact)
```

Explicit `to_html` and `to_svg` calls raise an `ImportError` with that install
hint when the renderer is absent. Notebook `_repr_html_` and `formatted()` do
not make Typst mandatory: they omit HTML and let the existing LaTeX or text
representation take over.

`TensorNetwork.__str__` remains its Graphviz DOT representation; use
`network.to_dot()` when that intent should be explicit. The rich-display
methods are separate and do not replace Symbolica's inherited `to_latex` API.

Expression and module-level HTML, SVG, and rich-display calls accept
`notation_source=...` for advanced customization. This is a trusted, complete
replacement for Tydenso's bundled `notation.typ`, not a style fragment; Typst
executes it and it must expose the notation interface expected by the bundled
renderer. Do not pass untrusted input. Display customization remains outside
the Atom payload, while portable representation and math-label declarations
continue to travel with the expression.

## Symbolic simplification

Community simplifiers operate on ordinary Symbolica expressions. Cross the
boundary explicitly, then reinfer when the result is still valid tensor syntax:

```python
from symbolica.community.idenso import simplify_metrics

simplified = simplify_metrics(dirac_trace.to_expression())
tensor_simplified = as_tensor(simplified)
```

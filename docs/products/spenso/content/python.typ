#import "../../shared.typ": boundary, callout, source-link, product-link

#let python = [
= Python tensor workflows

Spenso's Python interface is a native adapter over the same tensor structures, libraries, and
network executor as the Rust crates. Use the generated reference for exact signatures and this
guide for the object boundaries and execution sequence that those signatures do not explain.

== Compute a scalar product

Once the tensor module is available, this example computes the squared Euclidean
length of a vector. A representation defines its three component positions;
`Tensor.dense` supplies their values. Matching index labels request contraction.

// docs-example: compile
```python
from symbolica.community import tensor as sp

space = sp.Representation.euc(3)
p = sp.TensorName.vector("quickstart::p")(space)
vector = sp.Tensor.dense(p, [1.0, 2.0, 3.0])

network = vector("i") * vector("i")
result = network.to_tensor()
assert result.scalar() == 14.0
```

The product creates a computation; `to_tensor()` executes a copy and returns its
component data. Use `result[:]` for the flat components of a tensor with free
indices, or `result.scalar()` when no free indices remain. Call `network.execute()`
when you want to update the network itself and inspect its progress.

For symbolic calculations, work with `p` before supplying component values.
Multiplication automatically uses dot notation for two vector heads declared with
`TensorName.vector()`. A generic tensor or composite expression with one free
axis keeps its indexed contraction instead; `contract()` packs compatible vector
contractions into Schoonschip notation. Choose `dot()` explicitly or use
`to_dots()` to convert a compact scalar metric contraction.
For reusable numerical definitions, register `vector` in a `TensorLibrary` and
pass that library when evaluating an expression. Each method's reference entry
includes a runnable example followed by its parameter descriptions. The same
documentation is available in Python with `help(sp.Tensor.dense)`.

== Availability and version boundary

#boundary("A Symbolica community module", [
  Import Spenso as `symbolica.community.tensor`; there is no standalone `spenso` wheel. The
  published Symbolica wheel bundles this community module. Its Symbolica version determines
  which Spenso API it contains. A source checkout or generated `.pyi` file does not add the
  native module to an existing Python environment.
])

Install the current assembly and check the environment before building a workflow:

// docs-example: syntax
```sh
python -m pip install --upgrade symbolica
python -c "import symbolica.community.tensor as spenso; print(spenso.__name__)"
```

Source embedders can build a custom
#link("https://github.com/symbolica-dev/symbolica-community")[community-module assembly]. Record
the Symbolica assembly version with reproducible results; it is a more useful Python
compatibility fact than the version of an unrelated local Rust checkout.

Spenso, Idenso, and Symbolica core share one native library, one Symbolica kernel, and one
Python `Expression` type. The community assembly registers `SpensoModule` as
`symbolica.community.tensor_native`; its public wrapper imports the native exports and calls
the host-provided `initialize_module`. That initializer must remain in the export list.
`SpensoModule::get_name()` returns `tensor`; package the public wrapper under
`python/symbolica/community/tensor/`. Python classes and generated stubs use
`symbolica.community.tensor`, while Rust crate names and symbolic `spenso::` names remain unchanged.
Include the generated `__init__.pyi` beside that wrapper and a `py.typed` marker in
the `tensor` package so editors and type checkers can discover the public API.
The Spenso documentation exporter writes and checks the same stub in
`docs/api/python/spynso3.pyi` and the notebook host's `tensor` package.
Do not link Spynso into `gammaloop._gammaloop` or distribute it as a second native extension.

GammaLoop owns the Spynso source and bundled Tydenso `render.typ` and `notation.typ` assets;
Symbolica Community owns the wheel. Their native dependencies must resolve one Symbolica
source revision. Typst, its fonts, and offline packages are embedded in the Rust
extension. No Python compiler package is required.

=== Pyodide builds

The adapter's default `native` feature includes native arithmetic and C++ evaluator
compilation. For Pyodide, the community assembly must depend on `spynso3` with
`default-features = false` and forward its `wasm` feature to `spynso3/wasm`.
Every other Symbolica consumer in that assembly must also disable native defaults
and select its portable arithmetic features; Cargo combines features across dependencies.

`TensorEvaluator` remains available in this configuration. `TensorEvaluator.compile`
and `CompiledTensorEvaluator` require the `native` feature and are absent in Pyodide.

With Python 3.14, `pyodide-build==0.39.0`, `maturin==1.15.0`, and the `314.0.7`
cross-build environment installed, run the following from the community assembly's
directory. Its manifest must define the `release-small` profile and forward `wasm`
as described above; the GammaLoop workspace root builds a different Python package.

// docs-example: syntax
```sh
export RUSTUP_TOOLCHAIN="$(pyodide config get rust_toolchain)"
rustup target add --toolchain "$RUSTUP_TOOLCHAIN" wasm32-unknown-emscripten
pyodide build . --outdir dist --no-isolation \
  -C maturin.build-args="--locked --profile release-small \
  --no-default-features --features wasm \
  -- -C link-arg=-sEXPORTED_FUNCTIONS=_PyInit_core"
```

In this repository, `just check-spynso-wasm` checks the adapter for Emscripten and
rejects native arithmetic or workspace-hack dependencies in its portable graph.

== Choose the right object

Tensor factories, patterns, and name accessors use `dirac_gamma`, `color_f`, and
`color_t`: for example, `TensorExpression.dirac_gamma(4)`,
`TensorPattern.color_f(8, a_, b_, c_)`, and `TensorName.color_t()`.
The registered Symbolica heads remain `spenso::gamma`, `spenso::f`, and
`spenso::t`, including in saved expressions and printed mathematics.

The tensor package owns the full tensor vocabulary formerly exposed through
`hep.Symbols`. Use `TensorExpression.gamma0(d)`, `.charge_conjugation(d)`, and
`.levi_civita(rep, rank=4)` for open tensors with distinct ports; corresponding
`TensorPattern` and `TensorName` accessors support rewriting and head inspection.
Identical arguments to an antisymmetric head still vanish. The named expression
factories express distinct unresolved ports without changing that normalization.

Use `TensorPattern.dot/chain/trace` for compact syntax with wildcard operands or
factor sequences, and `PortPattern.chain_in/chain_out` for forward or reversed
contextual matrix channels. `TensorPattern.casimir/dynkin_index` also accept
arbitrary representation patterns. `Nc() -> Expression` constructs the canonical real
colour constant on demand, including its default numerical value 3. Importing the native
module therefore leaves time to call Symbolica's `set_license_key` before symbolic work.

`FactorProjector.symmetric`, `.antisymmetric`, and `.cyclic` group compatible
matrix factors for `chain` or `trace`. For example,
`chain(left, right, A, FactorProjector.antisymmetric(B, C), D)` projects just
the middle pair. Nested projector groups are supported; pass a precomposed
chain's individual factors separately. Symbolic groups stay compact and yield
`TensorExpression`; groups containing component tensors retain their data and
yield `TensorNetwork`. `expr.expand_projectors()` explicitly expands normalized
factor permutations while retaining the tensor interface and unrelated
factorization. The weights are `1/n!` for symmetric/antisymmetric permutations
and `1/n` for cyclic rotations.

- `Representation` and `Slot` define dimensions, duality, and abstract indices.
- `TensorExpression` carries a symbolic expression and its ordered tensor interface. A
  `Slot` supplies an explicit index; a `Representation` leaves a port open for later indexing.
- `TensorPattern` and `PortPattern` build rewrite patterns without requiring a concrete interface.
- `Tensor` owns dense or sparse data with a `TensorExpression` as its exact structure.
  Register a named tensor in a `TensorLibrary` for reuse by symbolic networks.
- `TensorNetwork` owns an expression graph. `ExecutionMode` selects the rewrite strategy:
  one smallest-degree rewrite per step, scalar work only, or the general smallest-degree strategy.
  Use `n_steps`, not the mode, to bound how many execution steps run.
- `TensorEvaluator` substitutes repeated numerical parameter batches into a symbolic tensor;
  its compiled counterpart owns the generated native evaluator.

The generated #link("reference/python/spynso3/Representation/")[`Representation`],
#link("reference/python/spynso3/Tensor/")[`Tensor`], and
#link("reference/python/spynso3/TensorNetwork/")[`TensorNetwork`] entries are the exact
versioned contracts. Do not infer index compatibility from two equal dimensions: names,
representations, and duality remain part of the value.

== Executable API tour and displays

#source-link("examples/notebooks/spenso_api_tour.py", label: "The complete Python API tour")
is a marimo notebook with one executed example for every declared type. Each entry shows its
construction, live display, Python representation, and declared members. Run it in a native
community-module environment with marimo and a C++ compiler on `PATH`:

// docs-example: syntax
```sh
python -m marimo edit examples/notebooks/spenso_api_tour.py
```

#source-link("examples/notebooks/spenso_notation.py", label: "The tensor notation showcase")
explains vectors and their family labels, free indices, unresolved and supplied ports,
scalar products, Dirac slashes, ordered chains, scoped copies, scalar invariants, and
component coordinates. It compares display settings and checks that compact contractions
preserve every component.

Mathematical objects expose notebook HTML and LaTeX displays. `TensorNetwork` draws its
current executable graph in notebooks, using Linnest's operator and typed-leaf styles.
The renderer passes native node, edge, and half-edge identities directly to Linnest's
graph builder; DOT remains a separate export format.
`render(config=...)` returns a displayable `DiagramRender` snapshot, while
`to_linnest(config=...)` returns a self-contained Typst document embedding the SVG.
Both accept typed `tensor.RenderSettings` values. Rust draws the graph directly
and the embedded Typst compiler typesets labels without graph plugins or MiTeX.

Use `drawing.to_svg()` or `drawing.to_html()` for string exports. Graph settings are
immutable typed values with discoverable constructors and properties:

```python
from symbolica.community.tensor import RenderSettings, LayoutSettings, StrokeStyle

drawing = network.render(config=RenderSettings(
    layout=LayoutSettings(layout_algo="dot"),
    edge_stroke=StrokeStyle(paint="#6f4d85", thickness=1.5),
))
drawing  # Displays directly in IPython, Jupyter, and Marimo.
```

`LayoutSettings`, `StrokeStyle`, and `DiagramRender` are also available from
`symbolica.community.hepkit`; both modules expose the same shared types.

`to_html(config=...)` wraps the graph in a figure labelled `TensorNetwork`.
Its execution summary uses `network.status` to show remaining nodes, operations,
contractions, and ready operations. “Graph reduced” means no graph work remains;
library or deferred results can still require materialization. Displaying the network
does not execute it; redisplay a stepped or executed network to update the snapshot.
Tensor labels use the registered tensor-name printer, including custom Typst heads
and scalar arguments. Edge slots use the same Typst index notation as the slot display, with
alphabet aliases assigned across the entire network. The inspector retains exact indices.
Hover summaries identify the node or edge kind. Clicking
shows tensor storage and ports, operator inputs, or contraction indices and metric
information; keyboard focus and Enter/Space use the same inspector.
`expression()` retains the semantic source formula, `to_dot()` exports the operation
graph, and `result_tensor()` shows evaluated component data. Rendering a network does
not execute it. Tensor expressions use ordinary multiplication when explicit indices
identify the slots, or when canonical multiplication preserves unresolved factor
occurrences in positional order. Brackets remain where sorting or combining factors
would lose that information, as in `p.outer(p)` with unresolved slots. Assigning explicit
indices removes brackets that are no longer needed. Nested product brackets flatten;
sums are not distributed and chain/trace scopes remain intact. Repeated indexed products
therefore become one n-ary network node before execution.
Settings and filters
print their full constructor arguments, and libraries and evaluators show concise summaries.
The policy types expose named constants and integer conversion; they are PyO3 classes, not
Python `enum.Enum` subclasses with `.name` and `.value` attributes.

The generated stubs distinguish symbolic results from concrete networks and type component
access as `Expression | float | complex`. For variadic `chain` and `trace`, symbolic-only
arguments retain `TensorExpression`; a concrete first or second factor selects
`TensorNetwork`. General mixed argument sequences retain the union return type. This avoids
unresolved types from overlapping variadic-tuple overloads in current type checkers.
`TensorExpression` subclasses Symbolica's `Expression` and retains the shared symbolic
tensor's ordered interface. Ordinary Symbolica inspection, matching, conversion, and
evaluation methods are inherited. Tensor algebra and checked rewrite overrides return
`TensorExpression`; other inherited methods retain Symbolica's return types. Tensor-specific methods
include `dirac_gamma`, which constructs a Dirac tensor, and `__getitem__`, which translates
between logical flat indices and coordinates.

The installed API regression checks compare every declared class member with the runtime,
exercise all tour examples, round-trip settings representations, and check real, complex,
sparse, and compiled evaluator results. The separate static fixture checks inferred return
types for indexing, arithmetic, composition, simplification, and evaluation.

== Metrics and oriented identities

`TensorExpression.g(left, right)` accepts a `Representation` or `Slot` for either
port. Slots supply their indices directly; representations leave ports unresolved.
Omitting `right` leaves the second port unresolved in the first port's representation.
For example, `TensorExpression.g(r("i"), r("j"))` constructs the indexed metric,
and `TensorExpression.g(r("i"), r)` leaves one index to fill. Its ports follow the
specified logical order. The representations must identify the same space and
have exactly equal dimensions; either port may carry its dual orientation.
For a fundamental color identity, pair the fundamental space with its dual:

// docs-example: compile
```python
from symbolica import E
from symbolica.community.tensor import Representation, TensorExpression

fund = Representation.cof(3)
identity = TensorExpression.g(fund, fund.dual())
indexed_identity = TensorExpression.g(fund("i"), fund.dual()("j"))
assert identity("i", "i").contract().to_expression() == E("3")
```

The same constructor supports symbolic dimensions. Distinct symbolic dimensions
are rejected, even when no numerical values have been assigned.

== Generated symbolic components

A named tensor with concrete dimensions needs no manually supplied component array.
If the library has no entry for it, network parsing generates symbolic components.
`TensorExpression.components(library=None)` exposes those same values in flat logical
row-major order, independent of the abstract index labels used in a contraction.

// docs-example: compile
```python
from symbolica.community.tensor import Representation, TensorName

spinor = Representation.bis(4)
J = TensorName("J")(spinor)
parameters = J.components()
assert len(parameters) == 4
assert J("s").components() == parameters
```

The values are the actual component expressions used by the parser, such as
`J(cind(0))`; use them directly as evaluator parameters or substitution keys.
The accessor executes a temporary network through the normal parser and executor.
For compound expressions it returns the contracted result; for a scalar it returns
a one-element list. It neither changes the expression nor registers generated tensors.
All enumerated dimensions must be concrete.

A supplied `TensorLibrary` is used for parsing, execution, and result extraction.
Registered tensors contribute their stored values instead of generated symbols.
The default is the ordinary HEP library; `TensorLibrary.hep_lib_atom()` selects
its atom-valued variant. Component values retain the same Python types as
`TensorNetwork.result_tensor()`: expressions, floats, or complex numbers.

== Construct and inspect concrete data

This complete source creates a named rank-two tensor, verifies one component, and converts its
storage in place. The documentation harness compiles the Python source without importing the
native module.

// docs-example: compile
```python
from symbolica.community.tensor import Representation, Tensor, TensorName

rep = Representation.euc(2)
i = rep("i")
j = rep("j")
structure = TensorName("A")(i, j)
matrix = Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0])

assert len(matrix) == 4
assert matrix[0, 0] == 1.0
matrix = matrix.to_sparse()
assert matrix[1, 1] == 1.0
```

`Tensor.dense` requires row-major data whose length is the product of the structure dimensions.
`Tensor.sparse` instead needs the element type and starts empty. `to_dense()` and `to_sparse()`
return independent tensors with the requested storage representation. Assign the result to
keep it; the original tensor retains its storage and values. They do not change slots or
re-index the tensor. See the
#link("reference/python/spynso3/Tensor/#exports-tensor-dense-associatedfunction")[dense constructor] and
#link("reference/python/spynso3/Tensor/#exports-tensor-to-sparse-method")[conversion contract].

#callout("Diagnose structure before storage", [
  A constructor failure usually means the data length and dimensions disagree. An unexpected
  contraction or exterior product is instead an index/duality problem. Print the structure and
  slots before changing dense/sparse storage, because a storage conversion cannot repair a
  structural mismatch.
])

== Transforming expressions and axes

Use ordinary algebra directly on a `TensorExpression`: `expand`, `expand_num`, `factor`,
`collect`, `collect_symbol`, `collect_num`, `collect_factors`, `collect_by_coefficient`,
`collect_horner`, `together`, `cancel`, and `apart` return `TensorExpression` values.
They preserve logical port order, unresolved port identities, and explicitly assigned
metadata. A zero result retains its rank and slots. `expand_num()` distributes numerical
coefficients over sums; it does not request full polynomial expansion. `collect_factors()`
extracts factors common to terms in sums, including nested sums.

// docs-example: compile
```python
from symbolica import Expression, S
from symbolica.community.tensor import Representation, TensorName

x, y = S("x", "y")
A = TensorName("algebra::A")(Representation.euc(3))("i")
expression = 2 * (x + y) * A
expanded = expression.expand_num()
factored = expanded.collect_factors()
assert factored.structure.axes == expression.structure.axes
assert (factored - expression).expand().to_expression() == 0
assert isinstance(expression, Expression)
assert expression.replace(x, y).structure.axes == expression.structure.axes
assert expression.derivative(x) == 2 * A
# Coefficient extraction inherits Symbolica's ordinary Expression return type.
assert type(expanded.coefficient(x)) is Expression
```

These operations use Symbolica's algebra implementation. Trusted rearrangements reuse
the established interface; user-supplied collection callbacks require checking the result.
A callback that removes or changes external ports raises an error. Label unresolved axes
before a transformation if their occurrence identities cannot be tracked through the result.

`replace(TensorRule(...))` applies checked whole-tensor replacement rules. The same method
accepts Symbolica's `replace(pattern, rhs, ...)` form, including callbacks and traversal
options. `replace_multiple`, `map`, and `derivative` also return tensors and check that the
result preserves the ordered interface. Derivatives of scalar coefficients retain typed
zeros; formal derivatives of unknown tensor functions require component differentiation
with `Tensor.map_components` or a tensor-name derivative callback.

Inherited operations such as `coefficient` and `terms` return ordinary `Expression`
values, which carry no tensor interface metadata. Tensor expressions are accepted directly
by Symbolica functions and matching APIs. Use `to_expression()` when deliberately selecting
the base expression behavior, for example to make a rewrite that changes tensor rank.
To rebuild a tensor after a structure-preserving raw operation, pass
`structure=original.structure` to the `TensorExpression` constructor; the constructor checks
the declared interface.

All three tensor forms expose `structure`, `axes`, `rank`, `shape`, and `is_scalar`.
The immutable `structure` is the canonical signature of the free axes; `axes` and
`shape` follow the current component or port order.
`index(...)` (also `obj(...)`) fills only unresolved ports. `reindex(...)` assigns every
external port and permits intentional contractions; `AUTO` leaves a port unchanged.
`rename_indices({...})` performs simultaneous substitutions on external indices and
rejects contractions or capture of internal dummy indices. Use immutable `Slot` keys
when the same label occurs in multiple representations. Tensor and network operations
keep their actual stored values; they do not resolve a fresh copy from a library.

`permute_axes([2, 0, 1])` reorders the public axes. On a concrete tensor it transposes
the data; on an expression or network it changes the external interface consistently.
Unresolved expression/network ports acquire fresh labels to keep their occurrences
distinct. Use `reindex(...)` afterwards to choose displayed labels.

Nested index payloads use `intern="indices"` for reversible encoding or
`intern="flattened"` for readable flattened names. Both modes affect only indices,
leaving tensor heads, scalar arguments, and dimensions intact. The default `intern=None`
leaves indices unchanged. Constructors, calls, `index`, `reindex`, and `rename_indices`
accept the same literal choices.

== Component operations

`tensor.dtype` is the actual component class: `float`, `complex`, or Symbolica's
`Expression`. `tensor.storage` is `"dense"` or `"sparse"`. `copy()` and Python's
`copy.copy()` duplicate concrete storage or network state.

Integer access uses flat logical row-major order, including negative indices.
Coordinate tuples follow `axes`. A flat slice returns a flat list;
coordinate slices such as `tensor[:, -1]` or `tensor[::-1, :]` return nested lists
for the selected axes. Assignment accepts an integer or full integer coordinates,
including negative indices. Slice assignment and ellipsis indexing are not supported.

`map_components(callback, dtype=None)` returns a new tensor. By default the output
keeps the input component type; specify `dtype=Expression`, `float`, or `complex`
to convert it. A sparse tensor's implicit zero is mapped once. If that result is
nonzero, storage becomes dense so subsequent contractions include those entries.
Callbacks should depend only on the component value; traversal order is unspecified.

`Tensor.from_numpy(expression, array)` checks the complete logical shape and copies
real numeric data into float64 storage or complex data into complex128 storage.
`tensor.to_numpy()` returns an independent array in logical axis order, including
for noncontiguous input arrays or permuted tensors. NumPy is imported only when these
methods are called. Symbolic components require evaluation before NumPy export.

== Execution and simplification

`expression.to_tensor(library=...)` parses, executes, and returns concrete component
storage in one call. `network.to_tensor(...)` executes a copy, preserving the source
network's current progress. Both accept a `function_library` for broadcast callbacks.
Use the existing mutating `execute(...)` and `result_tensor(...)` methods when progress
should remain in the network.

`network.status` returns an immutable `ExecutionStatus`: remaining nodes, operation
nodes, internal contraction edges, currently ready operation labels, and `complete`.
These are graph counts, not cost estimates. `network.step(...)` returns a copy after
one native Single step, including the executor's normal preprocessing. Self traces
can therefore run before a reported ready operation. Repeat `step`, inspect or display
the returned network, or execute it to completion.

`expression.contract(...)` performs structural contractions and optional
chain or trace collection without evaluating algebraic identities. Independently,
`expression.simplify_algebra(...)` selects gamma, color, and epsilon families
with their existing dimension and gamma5 conventions. Python enables gamma and
colour by default, while epsilon is opt-in;
`simplify_algebra(gamma=True, color=True, epsilon=True)` selects all three. Each identity includes its necessary local
contractions. No nested contraction settings or implicit global expansion is involved.
The default `collect_coefficients=True` combines generated coefficients in
contracted Dirac traces and their bound vector factors, revealing cancellations
and applying `expand_num()` numerical distribution inside those generated
coefficients while preserving unrelated input sums and scalar prefactors. Set it to `False`
to retain nested identity output; free traces remain factorized in either mode.
Other identity kernels keep their existing coefficient normalization.
The result is a `TensorExpression` and reports `Complete`, `Deferred`, or `Capped` through
`reduction_status`. A `max_steps_per_domain` budget caps work while retaining the exact
remaining expression. Explicit `expand()` requests polynomial materialization.
The default budget is `None` (unlimited). The
#product-link("idenso", page: "reference/interfaces/", label: "Idenso reduction guide")
explains contraction modes, distribution permissions, representation filters and
family-specific output settings with typed examples.

== Numerical evaluation

Interpreted and compiled tensor evaluators share `parameters`, `input_size`,
`output_shape`, `supports_real`, `evaluate`, and `evaluate_complex`. Parameters are
reported in the original order. Both evaluators reject any batch row with the wrong
length before entering a native evaluation kernel. An empty batch produces no tensors;
a constant tensor takes one empty row per requested evaluation.

Compilation generates a complex kernel and, when the exact coefficients permit it,
a separate real kernel. Real evaluation rejects complex coefficients; use
`evaluate_complex` to retain them. Every result keeps the original logical axes and
tensor identity. Compilation creates the requested source/library files plus real
companions when applicable.

== Typed tensor factories and patterns

Calling a user-defined `TensorName` places scalar key arguments before structural ports.
Predefined tensors instead have typed factories on `TensorExpression`: `g`, `flat`, `gamma`,
`gamma5`, `projm`, `projp`, `sigma`, `f`, and `t`. Factory dimensions select representations;
they are not scalar arguments and do not add fields to tensor-library keys.

// docs-example: compile
```python
from symbolica.community.tensor import _, TensorExpression

gamma = TensorExpression.dirac_gamma(4)
gamma_ijmu = gamma("i", "j", "mu")
line = gamma(_, _, "mu") * gamma(_, _, "nu")
dirac_trace = line.trace()

generator = TensorExpression.color_t(8, 3)
T_aij = generator("a", "i", "j")
```

Gamma's public argument order is its stored interface: bispinor-in, bispinor-out, then
Minkowski. `_` (also exported as `AUTO`) leaves a local port unresolved; it is not a shared
Einstein index. Calling a partially indexed expression assigns only its remaining open ports.
Unresolved axes display as hollow squares; several unresolved axes on one tensor carry
zero-based axis numbers. Scalar parameters do not count as axes, and matching square
numbers on different tensors do not imply contraction. Hover over a square in HTML output
to see its axis and representation. Assigned indices retain their usual alphabet notation.
Raw predefined `TensorName` accessors expose heads for inspection and matching, not concrete
construction.

Patterns retain ordinary Symbolica wildcard behavior, including wildcard dimensions:

// docs-example: compile
```python
from symbolica import S
from symbolica.community.tensor import TensorPattern

D_, i_, j_, mu_ = S("D_", "i_", "j_", "mu_")
gamma_pattern = TensorPattern.dirac_gamma(D_, i_, j_, mu_)
```

`TensorPattern` shortcuts follow the same index order as concrete factories. General patterns
place scalar `args` before structural `ports`; `PortPattern.exact`, `.any`, `.self_dual`, and
`.dualizable` select fixed or constrained representation heads. Wildcard-head names are
strings ending in one underscore; dimensions, indices, and scalar arguments are numbers or
Symbolica expressions. Patterns can be passed directly to Symbolica replacement operations.

== Dirac adjoints and external-leg labels

`TensorExpression.dirac_adjoint()` uses Idenso's conjugation rules and the
registered gamma-zero factors at open bispinor ports. Its default matrix
convention exchanges the two endpoint labels of each open chain. When squaring
an amplitude whose diagrams pair external fermions differently, use
`dirac_adjoint(preserve_indices=True)` to keep every label attached to the same
physical leg. This convention distributes over the whole amplitude, so callers
do not need to find and exchange endpoints separately in each channel.

// docs-example: compile
```python
from symbolica import E
from symbolica.community.tensor import TensorExpression

gamma = TensorExpression.dirac_gamma(4)
direct = gamma("i", "j", "mu").to_expression() * gamma("k", "l", "mu").to_expression()
crossed = gamma("i", "l", "nu").to_expression() * gamma("k", "j", "nu").to_expression()
amplitude = TensorExpression((1 + 2 * E("1𝑖")) * direct + (3 - E("1𝑖")) * crossed)
adjoint = amplitude.dirac_adjoint(preserve_indices=True)
restored = adjoint.dirac_adjoint(preserve_indices=True).to_expression()
assert (restored - amplitude.to_expression()).expand() == E("0")
```

Both conventions conjugate scalar coefficients and retain the required boundary
factors. The physical-leg convention distributes those factors across sums
before canceling them. The default keeps its existing factored output. To form
a squared amplitude, distinguish the adjoint's indices with `wrap_indices` and
supply the appropriate shared completeness tensors. `tensor.wrap_indices(S("bra"))`
returns a `TensorExpression` with native scoped indices, preserving both its
logical interface and internal contractions. The original and scoped indices
remain distinct until a completeness tensor connects them. Alphabet display
uses the same base letters with primes for the scoped copy; raw display retains
`spenso::index_scope(scope, index)`. Repeating the same outer scope is idempotent,
and different scopes can nest. No index cooking is required. Keeping labels attached to
legs does not itself perform a spin sum. The new keyword requires a community
assembly built with this Spynso version; changing a stub file does not update
the native extension.

Gamma simplification also respects matrix orientation. In a compact chain,
`gamma(in,out,...)` is an ordinary matrix and `gamma(out,in,...)` is transposed.
Reversed gamma, gamma-zero, gamma-five and chiral-projector factors remain
explicit instead of entering forward-matrix Clifford or projector rules.
Ordinary subwords can still simplify, but mixed ordinary/transposed words are
not treated as if all their factors had the same orientation. Canonical chain
ordering does not supply a general transpose algebra. A supported explicit
four-dimensional `spenso::charge_conjugation` sandwich is an exception: Idenso
transposes each known gamma, gamma-zero, gamma-five or projector factor with
the sign fixed by $C = -i γ^2 γ^0$. Factor order is preserved and scalar
coefficients are not conjugated. Gamma/slash arguments of symbolic dimension
remain opaque. The same chain collector joins common starts or ends only in
self-dual spaces, reversing one word and exchanging its `in`/`out` markers.

Numerical evaluation with the shared HEP matrix data remains a separate way to
check such expressions. An installed-host comparison validates all sixteen
matrix components of a mixed ordinary/transposed product against that data.
The new `installed_charge_conjugation.py` regression additionally compares C
sandwiches with independent products of those Weyl matrices. Its initial
symbolic identities pass in the rebuilt host; full installed-host validation
awaits fixture setup corrections.

== Register data and execute a network

A symbolic network resolves named tensors through a library. Register the data first, construct
the expression with the same name and compatible slots, then execute and request the result:

// docs-example: compile
```python
from symbolica.community.tensor import (
    ExecutionMode,
    Representation,
    Tensor,
    TensorLibrary,
    TensorName,
    TensorNetwork,
)

rep = Representation.euc(2)
A = TensorName("A")
structure = A(rep, rep)
library = TensorLibrary()
library.register(Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0]))

network = TensorNetwork(A(rep("i"), rep("j")), library=library)
network.execute(library=library, mode=ExecutionMode.All)
result = network.result_tensor(library=library)
assert len(result) == 4
```

The library also acts as a mapping from full unresolved signatures to stored
component data. `library[structure]` returns a `Tensor`; its `expression()` method
returns the symbolic reference. Keys include scalar arguments and ordered
representations, so tensors such as `A(x, 7, rep)` and `A(x, 8, rep)` can coexist.
`library["A"]` or `library[A]` (with `A` a `TensorName`) is a convenience when the
name identifies exactly one stored signature; ambiguous or absent names raise
`KeyError`.

// docs-example: compile
```python
stored = library[structure]
reference = stored.expression()
signatures = library.keys()
assert structure in library
assert library.get("missing_tensor") is None
for signature, tensor in library.items():
    print(signature, tensor.structure.shape)
```

`keys()`, `values()` and `items()` return snapshot lists in matching order;
iteration yields signatures and `len(library)` counts stored tensors. Returned
tensors are independent copies: register an edited tensor again to replace the
stored data. Dimension-dependent factories, such as metrics, can be accessed
with an exact concrete signature but do not appear among stored entries.

`key in library` checks whether lookup can resolve the signature, including
dimension-dependent factories, without constructing component data.
`library.get(key, default=None)` returns an independent tensor or the supplied
default when the signature is absent. Both accept the same keys as indexing;
ambiguous names and invalid keys still raise errors. In particular, a factory
signature can be in an otherwise empty library while `len(library)` remains zero.

Displaying a library in a notebook opens a compact catalogue of mathematical
signatures. Select a tensor to inspect its components in the usual Memory grid
or Matrix view. The collapsed Python panel provides executable access and
expression examples, including full namespaces, scalar arguments and custom
representations. Its Print tab shows text, Typst and LaTeX output calls.
`library.to_html()` returns the same standalone HTML; pass `settings` to use
the existing display options. Large catalogues show a bounded preview while
`keys()` and `items()` continue to expose every stored entry.

`ExecutionMode.Single` selects the smallest-degree single-rewrite strategy, but execution still
continues while work remains unless it is bounded; use `n_steps=1` to inspect exactly one step.
`Scalar` processes scalar work while retaining tensor structure, and `All` attempts the complete
available execution. A successful `execute()` does not make `result_scalar()` appropriate for a tensor result; choose
`result_tensor()` or `result_scalar()` according to the remaining structure. Exact behavior is
linked from #link("reference/python/spynso3/TensorNetwork/#exports-tensornetwork-execute-method")[`execute`],
#link("reference/python/spynso3/TensorNetwork/#exports-tensornetwork-result-tensor-method")[`result_tensor`],
and #link("reference/python/spynso3/ExecutionMode/")[`ExecutionMode`].

#callout("Interpret a network failure by ownership", [
  A missing-library error means a symbolic tensor name has no registered data. A result-kind
  error means the graph still has tensor structure when a scalar was requested, or vice versa.
  A changed network structure invalidates assumptions about contraction order; replacing only
  values with the same registered structure does not.
])

== Mathematical display

Concrete tensors open an interactive component explorer by default in notebooks.
Choose two displayed axes, fix the other coordinates, and compare slices.
The visible *Memory grid / Matrix* toggle switches the current display for either
one slice or all slices. Axis and slice controls start collapsed into a settings
summary; expand the chevron to edit them. The soft fields share one row when
space permits and wrap into two columns in narrow notebook outputs.
Memory-grid shading always reflects component bytes. Selecting a component
shows its formula through the same tensor printer, including custom names.
The explorer follows the notebook's light or dark theme and fills its output
width. Matrix formulas scale to fit their cells; select a cell to inspect the
formula at full size above the matrix.
Deeper blue indicates a larger component payload: encoded Symbolica expression
bytes, or fixed-width numeric storage. Sparse implicit entries share one default;
they do not each occupy a stored component. Counts exclude allocation overhead
and shared symbol metadata. Large tensors carry a bounded preview without
densifying sparse storage; unloaded cells are explicitly marked rather than
shown as zero. Heaviest-first ordering applies to the included components.

For a static notebook matrix, use
`tensor.formatted(settings=DisplaySettings(tensor_view="matrix"))`.
The same setting works with `to_html`. SVG, Typst and LaTeX retain their
mathematical matrix/slice output. Symbolic tensor expressions retain their
existing display and index alphabets. The explorer runs inside a self-contained
sandboxed frame and requires no live Python callbacks or external assets.

`TensorExpression` and `Tensor` expose semantic display methods. A network's
`expression()` provides the same formula display separately from its graph.
Symbolic notebook outputs render the complete expression in HTML, LaTeX, and
pretty text, regardless of its size. Display preserves factorization and ordered
ports without expanding or simplifying the calculation.

The notebook frontend can impose its own output limit. Marimo's
`tool.marimo.runtime.output_max_bytes` setting in `pyproject.toml` controls that
separate limit; a frontend warning does not mean the tensor renderer omitted
terms.

`TensorExpression.formatted()` and the free `formatted(expression)` defer each
backend until the notebook requests it and cache that representation. Their
plain-text fallback also includes the complete expression. Settings are captured
when creating the display value; custom printers run when a backend is first
requested. Create a new formatted value to reflect later changes to printer rules.
The explicit `to_html()`, `to_svg()`, `to_latex()`, `to_typst()`, and
`format_tensor()` methods provide the same complete expression in the requested
format.

`DisplaySettings` controls the ports, Schoonschip, and call layouts, dimensions, parentheses,
commas, symbol scripts, component notation, and index/factor spacing. Positional calls such as `to_typst(True)`
and `formatted(True)` still request dimensions. Rich Typst output collects inverse factors
at the same product level into a single fraction, including rational coefficients.
The default `ports` layout displays vectors inserted into tensors with bras and kets.
Keep the default settings, or select `DisplaySettings.ports()`, for this notation.
`DisplaySettings.schoonschip()` selects the alternative compact layout with bold momentum
labels in the tensor's index positions. Both layouts retain unresolved AUTO ports.
Schoonschip notation writes a vector in the position of an index contracted with it:
$T(p, nu) = T_(mu nu) p^mu$. See
#link("https://www.nikhef.nl/~form/maindir/documentation/tutorial/book.pdf#page=14")[A. Heck,
_FORM for Pedestrians_, §1.2.2, pp. 9–10]. For the metric this gives
$g(p, q) = g_(mu nu) p^mu q^nu = p dot q$. The notation showcase renders the indexed
metric and vectors, the default bra-and-marker display, and the dot product directly
from Spenso objects. Its foldout shows the generated Typst source and alternative
layouts with vector arguments or vectors in index positions.
`to_dots()` converts the compact metric form without evaluating component data.
Compound graph indices and generated dummy indices use each representation's alphabet by default: Lorentz indices
render as $mu, nu, rho, sigma$, fundamental color as $i, j, k, l$, and bispinor or adjoint
color as $a, b, c, d$. The alphabet repeats with subscripts when needed. Repeated indices
share a label throughout an expression, and existing numeric or manually named labels are
reserved to prevent collisions. Scoped copies keep the same base letters with primes;
internal dummy identifiers never become alphabet labels.

Color invariants use compact notation by default: $C_F$, $C_A$, $T_R$, and $N_c$.
Select `DisplaySettings(invariant_style="explicit")` to display the degree and
representation as $C_2(F)$, $C_2(A)$, and $I_2(F)$. Higher invariants always use
$C_k(R)$, $I_k(R)$, and $G_k(R,S)$, with $F$ and $A$ denoting the fundamental and
adjoint representations. These scalar arguments retain parentheses and commas
independently of the tensor-index layout settings.

`show_dimensions=True` adds representation dimensions and disables compact
quadratic aliases, for example $C_2(F_3)$ and $C_2(A_8)$. Use it to distinguish
multiple color spaces. Older scalar constants carry no dimension metadata.
These settings affect presentation only: exact expressions and representation
identities are preserved across plain text, LaTeX, Typst source, HTML and SVG.

// docs-example: compile
```python
from symbolica.community.tensor import DisplaySettings, Representation, TensorExpression

invariant = TensorExpression(Representation.cof(3).casimir())
explicit = DisplaySettings(invariant_style="explicit")
source = invariant.to_typst(settings=explicit)  # C_2(F)
latex = invariant.to_latex(settings=explicit)
preview = invariant.formatted(settings=explicit)
```

Graph-derived indices record their origin as an edge, half-edge or vertex identifier
and a local index label. Choose `DisplaySettings(index_style="graph")` to display this
origin in abbreviated form, for example
$mu_(upright("h4"))$ for `hedge(4,1)`, $mu_(upright("e4"))$ for `edge(4,1)`, or $mu_(upright("v4"))$ for `vertex(4,1)`.
The local index label is implicit when it equals one; other values remain visible, as in
$mu_(upright("h4.2"))$ for `hedge(4,2)`.
Use `index_style="alphabet"` for the compact default or `index_style="raw"` for the original
symbolic index notation. Changing the display style preserves the underlying expressions, tensor
interfaces, and exact notebook payloads. In raw mode, endpoint labels retain their existing
subscript notation, including distinct higher-spin and dummy slots.
After `to_expression()`, ordinary Symbolica printing owns namespace elision and nested
bracket highlighting. Use `format(show_namespaces=True)` to display qualified names.

Tensor names accept backend-specific head formatting through `print`:

// docs-example: compile
```python
from symbolica.community.tensor import Representation, TensorName

Jbar = TensorName(
    "Jbar",
    print={"typst": "macron(J)", "latex": r"\bar{J}"},
)(Representation.bis(4))
```

Mapping values are trusted backend source without math delimiters. Spenso adds ordinary
arguments, abstract indices, and component coordinates, including inside concrete tensor
matrices. An omitted backend keeps the ordinary name; `plain` can also be supplied.
A callable retains Symbolica's `print(expression, mode=..., **options)` signature and
replaces the complete tensor display. Returning `None` selects standard tensor notation.
Notebook rendering evaluates local callbacks before sending their visual output to Typst;
the callable itself is not serialized into portable Atom payloads.

Concrete tensor components default to $A(x,7)^(0,1)$. Choose
`DisplaySettings(component_style="array")` for $A(x,7)[0,1]$, or
`component_style="superscript"` for the default. Both styles keep ordinary
arguments in parentheses and apply to individual components and matrix entries.
Use `component.formatted(settings=DisplaySettings(component_style="array"))`
to show array notation in a notebook. The setting changes presentation only;
component coordinates remain independent of abstract-index styles.

// docs-example: compile
```python
from symbolica.community.tensor import DisplaySettings, TensorExpression

trace = TensorExpression.gamma5(4).trace()
source = trace.to_typst(
    settings=DisplaySettings(show_dimensions=True, parentheses=False)
)
```

// docs-example: compile
```python
compact_indices = numerator.formatted(settings=DisplaySettings(index_style="alphabet"))
graph_indices = numerator.formatted(settings=DisplaySettings(index_style="graph"))
raw_indices = numerator.formatted(settings=DisplaySettings(index_style="raw"))
```

`to_typst` and `format_tensor` emit source using the ports layout. Schoonschip, call, and
custom-spacing settings require Tydenso's Typst notation layer; use HTML, SVG, or notebook
display for those settings. Source-only methods reject unsupported settings rather than
silently ignoring them.

A supplied rank-one tensor fills an argument position in the ports layout.
A filled marker distinguishes it from an unresolved axis: inline-metric spaces
use `■`, self-dual spaces use `●`, and a dualizable space and its dual use
left- and right-pointing filled triangles. A hollow square `□`
marks an unresolved axis. The supplied tensors appear as bras or kets according
to their representation; this notation does not conjugate their components.

A momentum contracted into a gamma matrix displays as a slash, with its two
remaining bispinor positions. Gamma uses the same endpoint notation as generic
tensors: a supplied row spinor becomes a bra and occupies a filled position.
Explicit spinor indices remain visible; unresolved ones use hollow squares.
The consumed Lorentz position has no placeholder because the slash encodes it.
`contract()` collects connected matrix factors even after a rank-one tensor
has filled an endpoint. The supplied endpoint does not add an external axis:
for example, a spinor times two slashed gamma matrices retains one open spinor
slot. Collected chains retain matrix order and the same endpoint markers;
collecting a trace does not evaluate it. Set `collect_chains=False` to retain
the individual matrix factors.

// docs-example: compile
```python
from symbolica.community.tensor import Representation, TensorExpression, TensorName

mink = Representation.mink(4)
p = TensorName.vector("p")
indexed = TensorExpression.dirac_gamma(4)("a", "b", "mu") * p(mink("mu"))
pslash = indexed.contract()
pslash.formatted()
```

HTML and SVG rendering use the embedded compiler:

// docs-example: compile
```python
from symbolica.community.tensor import DisplaySettings, TensorExpression

trace = TensorExpression.gamma5(4).trace()
compact = DisplaySettings.schoonschip()
html = trace.to_html(settings=compact)
svg = trace.to_svg(settings=compact)
rich = trace.formatted(settings=compact)
```

Python uses the bundled Typst render/notation assets directly, without calling the Tydenso
Wasm plugin. Explicit `to_html` and `to_svg` calls report compilation errors.
Notebook `_repr_html_` and `formatted()` retain their LaTeX or text fallback on rendering errors. `TensorNetwork.__str__` prints the source formula and `to_dot()` returns
the current graph. Its graph renderer uses Linnet's prepared-render pipeline and reports
rendering errors directly. Use `to_expression().to_latex()` when raw Symbolica notation is
needed instead of tensor-aware notation.

HTML output keeps selectable native MathML and embeds the same STIX Two Math font as
the documentation site, including in standalone offline notebooks. Math defaults to
21 px with padding and horizontal scrolling for long expressions. Set the CSS custom
property `--spenso-math-font-size` on the notebook or an output container to change
the base size; nested scripts retain their relative sizing.

Idenso reductions return typed `TensorExpression` values with these display methods.
Call `to_expression()` explicitly when an ordinary Symbolica expression is needed.
HTML, SVG, and rich-display functions accept `notation_source` as a trusted, complete
replacement for the bundled `notation.typ`, not a style fragment. Typst executes that source;
it must implement the expected notation interface and must not come from untrusted input.
Display customization stays outside Atom payloads; portable representation and math-label
declarations continue to travel with the expressions.

== Citation tracking

After a calculation, `symbolica.get_citations()` returns citation objects with
references, explanations in `reasons`, and bibliography entries from
`citation.to_bibtex()`. Spenso is authored by Lucien Huber; Idenso is authored by
Lucien Huber and Ben Ruijl.

Successful tensor operations cite Spenso. Symbolic contraction, algebra
simplification, canonicalization, Dirac adjoints, and index/notation rewrites add
Idenso. Every Idenso citation also credits Spenso for the tensor-expression
structure and display that define tensors in Symbolica. Component-network
execution and numerical evaluation add their own reasons to the Spenso citation.

Dirac gamma algebra also cites
#link("https://arxiv.org/abs/1203.6543")[FORM] by J. Kuipers, T. Ueda,
J. A. M. Vermaseren, and J. Vollinga. Color algebra also cites
#link("https://arxiv.org/abs/hep-ph/9802376")[Group theory factors for Feynman diagrams]
by T. van Ritbergen, A. N. Schellekens, and J. A. M. Vermaseren, the paper
underlying FORM's `color.h` package.

Citations accumulate during a session. Repeated operations keep a single entry
for each contribution.

== Repeated symbolic evaluation

Use `Tensor.evaluator()` when the tensor structure stays fixed and only symbolic parameters
change across batches. Supply constants, custom functions, and parameters explicitly. Call
`evaluate()` for real inputs and `evaluate_complex()` when coefficients or inputs are complex.
`compile()` creates source and a shared library, so treat its filenames, compiler, architecture,
and optimization level as reproducibility inputs rather than incidental arguments.

Control Symbolica-backed Rayon work with `SymbolicParallelism`: `Serial` keeps work on the
calling thread, `Auto` checks license capability and applies workload heuristics, and `Parallel`
forces Rayon without that safety choice. Configure the policy before benchmarking; otherwise a
threading-policy change can be mistaken for an algorithmic improvement.

For exact parameters and defaults, use the
#link("reference/python/spynso3/Tensor/#exports-tensor-evaluator-method")[evaluator reference] and
#link("reference/python/spynso3/set_symbolica_rayon_enabled-function/")[symbolic parallelism reference].
The implementation starts in
#source-link("crates/spynso3/src/lib.rs", label: "the Spenso Python adapter").
]

== Tensor evaluation and Symbolica API parity

`tensor.evaluator(params, functions=..., jit_compile=True, ...)` uses the same
argument order, defaults, function definitions and optimisation settings as
`Expression.evaluator`. It builds one Symbolica evaluator for all components in
logical order. `evaluate` and `evaluate_complex` accept the same array-like inputs
as Symbolica and return one `Tensor` per input row, preserving axes, name and
arguments. Sparse input and permuted axes follow the same component ordering.
SymJIT compilation occurs on first numerical use; `jit_compile=False` selects
interpreted evaluation.

The former `constants, funs, params` interface is removed. Substitute fixed values
before constructing the evaluator, and pass `FunctionDefinition` objects through
`functions`. `TensorEvaluator.compile` now follows `Evaluator.compile`, including
its required `number_type`, compiler flags and native/SIMD/CUDA options. The
compiled object's `evaluate` method uses the selected number type. The old
compiled `evaluate_complex` method is removed. `scalar_evaluator` exposes the
underlying Symbolica object for arbitrary-precision evaluation, instruction
inspection, export and other scalar operations; these return component arrays
without tensor metadata.

An audit against Symbolica 3.0.1 found matching signatures for the explicitly
wrapped `TensorExpression` algebra operations: `apart`, `cancel`, `collect`,
`collect_by_coefficient`, `collect_factors`, `collect_horner`, `collect_num`,
`collect_symbol`, `derivative`, `expand`, `expand_num`, `factor`, `map`,
`replace_multiple` and `together`. Their tensor results additionally validate the
external interface. `TensorExpression.evaluator` is inherited from `Expression`
and evaluates its symbolic payload; use `to_tensor().evaluator(...)` to evaluate
materialized components.

The remaining differences are:

- `TensorNetwork.replace` still uses `level_range` and optional booleans, lacks
  `partial`, `once`, `bottom_up` and `nested`, and rejects replacement callbacks.
  `TensorExpression.replace` already delegates those controls to Symbolica;
  making its `rhs` optional is intentional, to accept whole-tensor rules.
- `TensorNetwork.evaluate` still accepts real constants and a map of Python
  functions. `Expression.evaluate` accepts real or complex substitutions and an
  optional decimal precision. The network route has different capabilities and
  does not yet provide scalar API parity.
- Tensor indexing, powers and calls carry tensor semantics. Tensor display methods
  add representation and layout settings; `to_typst` takes `show_dimensions`
  where the scalar method takes `show_namespaces`. These are intentional
  specialisations, not interchangeable scalar operations.
- Inherited scalar inspection and transformation methods can return ordinary
  `Expression` values. Use the explicitly wrapped algebra methods when the result
  must retain a checked tensor interface.

The audit also reproduced a shared Symbolica 3.0.1 C++ export failure at kernel
revision `942bd2c` for a non-real coefficient: compiling `Symbol.I * x` with `number_type="complex"`
generates `T(T(0e0), T(1e0))`, which is invalid for `std::complex<double>`.
Both scalar and tensor evaluators fail, with either direct-translation setting.
The older tensor-only exporter avoided this kernel bug through its own complex
coefficient printer. Delegating native compilation exposes it; the existing
complex-coefficient compilation regression remains failing until the shared
exporter is fixed. Interpreted and SymJIT evaluation, including the current
notebook, work. This is a backend defect, not an argument-signature mismatch.

`crates/spynso3/tests/installed_evaluator_parity.py` checks the shared signatures
and compares component results for function definitions, JIT and interpreted
execution, sparse data, permuted axes and constant scalars. Native real/complex
compilation is covered by `installed_api_operations.py`.

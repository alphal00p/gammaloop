#import "../../shared.typ": boundary, callout, source-link

#let python = [
= Python tensor workflows

Spenso's Python interface is a native adapter over the same tensor structures, libraries, and
network executor as the Rust crates. Use the generated reference for exact signatures and this
guide for the object boundaries and execution sequence that those signatures do not explain.

== Availability and version boundary

#boundary("A Symbolica community module", [
  Import Spenso as `symbolica.community.spenso`; there is no standalone `spenso` wheel. The
  published Symbolica wheel bundles this community module. Its Symbolica version determines
  which Spenso API it contains. A source checkout or generated `.pyi` file does not add the
  native module to an existing Python environment.
])

Install the current assembly and check the environment before building a workflow:

// docs-example: syntax
```sh
python -m pip install --upgrade symbolica
python -c "import symbolica.community.spenso as spenso; print(spenso.__name__)"
```

Source embedders can build a custom
#link("https://github.com/symbolica-dev/symbolica-community")[community-module assembly]. Record
the Symbolica assembly version with reproducible results; it is a more useful Python
compatibility fact than the version of an unrelated local Rust checkout.

Spenso, Idenso, and Symbolica core share one native library, one Symbolica kernel, and one
Python `Expression` type. The community assembly registers `SpensoModule` as
`symbolica.community.spenso_native`; its public wrapper imports the native exports and calls
the host-provided `initialize_module`. That initializer must remain in the export list.
Do not link Spynso into `gammaloop._gammaloop` or distribute it as a second native extension.

GammaLoop owns the Spynso source and bundled Tydenso `render.typ` and `notation.typ` assets;
Symbolica Community owns the wheel. Their native dependencies must resolve one Symbolica
source revision. The `gammaloop[typst-display]` extra adds only the optional renderer, not
another Spynso binary.

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
community-module environment with marimo, Typst, and a C++ compiler on `PATH`:

// docs-example: syntax
```sh
python -m marimo edit examples/notebooks/spenso_api_tour.py
```

Mathematical objects expose notebook HTML and LaTeX displays. `TensorNetwork` draws its
current executable graph in notebooks, using Linnest's operator and typed-leaf styles.
The renderer passes native node, edge, and half-edge identities directly to Linnest's
graph builder; DOT remains a separate export format.
`render(config=...)` returns interactive SVG and `to_linnest(config=...)` returns its
Typst entrypoint; both accept `linnet.RenderConfig`, like Feynman diagrams.
`to_html(config=...)` wraps the graph in a figure labelled `TensorNetwork`.
`expression()` retains the semantic source formula, `to_dot()` exports the operation
graph, and `result_tensor()` shows evaluated component data. Rendering a network does
not execute it. Settings and filters
print their full constructor arguments, and libraries and evaluators show concise summaries.
The policy types expose named constants and integer conversion; they are PyO3 classes, not
Python `enum.Enum` subclasses with `.name` and `.value` attributes.

The generated stubs distinguish symbolic results from concrete networks and type component
access as `Expression | float | complex`. For variadic `chain` and `trace`, symbolic-only
arguments retain `TensorExpression`; a concrete first or second factor selects
`TensorNetwork`. General mixed argument sequences retain the union return type. This avoids
unresolved types from overlapping variadic-tuple overloads in current type checkers.
`TensorExpression` deliberately specializes a few
inherited Symbolica names: `gamma` constructs a Dirac tensor, and `__getitem__` translates
between logical flat indices and coordinates. Their narrow override annotations document
this specialization instead of suppressing diagnostics for the entire module.

The installed API regression checks compare every declared class member with the runtime,
exercise all tour examples, round-trip settings representations, and check real, complex,
sparse, and compiled evaluator results. The separate static fixture checks inferred return
types for indexing, arithmetic, composition, simplification, and evaluation.

== Metrics and oriented identities

`TensorExpression.g(left, right)` creates an unresolved metric with ports in the
specified logical order. The representations must identify the same space and
have exactly equal dimensions; either port may carry its dual orientation.
Omitting `right` uses `left` for both ports, as in a Minkowski metric.
For a fundamental color identity, pair the fundamental space with its dual:

// docs-example: compile
```python
from symbolica import E
from symbolica.community.spenso import Representation, TensorExpression

fund = Representation.cof(3)
identity = TensorExpression.g(fund, fund.dual())
indexed_identity = identity("i", "j")
assert identity("i", "i").simplify_metrics().to_expression() == E("3")
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
from symbolica.community.spenso import Representation, TensorName

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
from symbolica.community.spenso import Representation, Tensor, TensorName

rep = Representation.euc(2)
i = rep("i")
j = rep("j")
structure = TensorName("A")(i, j)
matrix = Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0])

assert len(matrix) == 4
assert matrix[0, 0] == 1.0
matrix.to_sparse()
assert matrix[1, 1] == 1.0
```

`Tensor.dense` requires row-major data whose length is the product of the structure dimensions.
`Tensor.sparse` instead needs the element type and starts empty. `to_dense()` and `to_sparse()`
mutate the storage representation; they do not change slots or re-index the tensor. See the
#link("reference/python/spynso3/Tensor/#exports-tensor-dense-associatedfunction")[dense constructor] and
#link("reference/python/spynso3/Tensor/#exports-tensor-to-sparse-method")[conversion contract].

#callout("Diagnose structure before storage", [
  A constructor failure usually means the data length and dimensions disagree. An unexpected
  contraction or exterior product is instead an index/duality problem. Print the structure and
  slots before changing dense/sparse storage, because a storage conversion cannot repair a
  structural mismatch.
])

== Typed tensor factories and patterns

Calling a user-defined `TensorName` places scalar key arguments before structural ports.
Predefined tensors instead have typed factories on `TensorExpression`: `g`, `flat`, `gamma`,
`gamma5`, `projm`, `projp`, `sigma`, `f`, and `t`. Factory dimensions select representations;
they are not scalar arguments and do not add fields to tensor-library keys.

// docs-example: compile
```python
from symbolica.community.spenso import _, TensorExpression

gamma = TensorExpression.gamma(4)
gamma_ijmu = gamma("i", "j", "mu")
line = gamma(_, _, "mu") * gamma(_, _, "nu")
dirac_trace = line.trace()

generator = TensorExpression.t(8, 3)
T_aij = generator("a", "i", "j")
```

Gamma's public argument order is its stored interface: bispinor-in, bispinor-out, then
Minkowski. `_` (also exported as `AUTO`) leaves a local port unresolved; it is not a shared
Einstein index. Calling a partially indexed expression assigns only its remaining open ports.
Raw predefined `TensorName` accessors expose heads for inspection and matching, not concrete
construction.

Patterns retain ordinary Symbolica wildcard behavior, including wildcard dimensions:

// docs-example: compile
```python
from symbolica import S
from symbolica.community.spenso import TensorPattern

D_, i_, j_, mu_ = S("D_", "i_", "j_", "mu_")
gamma_pattern = TensorPattern.gamma(D_, i_, j_, mu_)
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
from symbolica.community.spenso import TensorExpression

gamma = TensorExpression.gamma(4)
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
supply the appropriate shared completeness tensors. Keeping labels attached to
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
from symbolica.community.spenso import (
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
for signature, tensor in library.items():
    print(signature, tensor.structure.shape)
```

`keys()`, `values()` and `items()` return snapshot lists in matching order;
iteration yields signatures and `len(library)` counts stored tensors. Returned
tensors are independent copies: register an edited tensor again to replace the
stored data. Dimension-dependent factories, such as metrics, can be accessed
with an exact concrete signature but do not appear among stored entries.

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
`DisplaySettings` controls the ports, Schoonschip, and call layouts, dimensions, parentheses,
commas, symbol scripts, component notation, and index/factor spacing. Positional calls such as `to_typst(True)`
and `formatted(True)` still request dimensions. Rich Typst output collects inverse factors
at the same product level into a single fraction, including rational coefficients.
Compound graph indices use each representation's alphabet by default: Lorentz indices
render as $mu, nu, rho, sigma$, fundamental color as $i, j, k, l$, and bispinor or adjoint
color as $a, b, c, d$. The alphabet repeats with subscripts when needed. Repeated indices
share a label throughout an expression, and existing numeric or manually named labels are
reserved to prevent collisions.

Choose `DisplaySettings(index_style="graph")` to retain graph identifiers, for example
$mu_(upright("h4"))$ for `hedge(4,1)`, $mu_(upright("e4"))$ for `edge(4,1)`, or $mu_(upright("v4"))$ for `vertex(4,1)`.
The second index is implicit when it equals one; other values remain visible, as in
$mu_(upright("h4.2"))$ for `hedge(4,2)`.
Use `index_style="alphabet"` for the compact default or `index_style="raw"` for the original
symbolic index notation. Graph-index styles preserve the underlying expressions, tensor
interfaces, and exact notebook payloads. In raw mode, endpoint labels retain their existing
subscript notation, including distinct higher-spin and dummy slots.
After `to_expression()`, ordinary Symbolica printing owns namespace elision and nested
bracket highlighting. Use `format(show_namespaces=True)` to display qualified names.


Tensor names accept backend-specific head formatting through `print`:

// docs-example: compile
```python
from symbolica.community.spenso import Representation, TensorName

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
from symbolica.community.spenso import DisplaySettings, TensorExpression

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


Contract a momentum with a gamma matrix, then collect its bispinor channel to display
slash notation. Chain collection keeps the ordered factors and leaves Dirac traces
unevaluated. The result remains a `TensorExpression` with its two bispinor ports.

// docs-example: compile
```python
from symbolica.community.spenso import Representation, TensorExpression, TensorName

mink = Representation.mink(4)
p = TensorName.vector("p")
indexed = TensorExpression.gamma(4)("a", "b", "mu") * p(mink("mu"))
pslash = indexed.schoonschip_net().collect_gamma_chains()
pslash.formatted()
```

Install the optional compiler to render HTML and SVG:

// docs-example: syntax
```sh
pip install 'gammaloop[typst-display]'
```

// docs-example: compile
```python
from symbolica.community.spenso import DisplaySettings, TensorExpression

trace = TensorExpression.gamma5(4).trace()
compact = DisplaySettings.schoonschip()
html = trace.to_html(settings=compact)
svg = trace.to_svg(settings=compact)
rich = trace.formatted(settings=compact)
```

Python uses the bundled Typst render/notation assets directly, without calling the Tydenso
Wasm plugin. Explicit `to_html` and `to_svg` calls raise an install-guidance `ImportError`
when the compiler is absent. Notebook `_repr_html_` and `formatted()` fall back to existing
LaTeX or text. `TensorNetwork.__str__` prints the source formula and `to_dot()` returns
the current graph. Its graph renderer uses Linnet's prepared-render pipeline and reports
rendering errors directly. These display methods do not replace Symbolica's inherited
`to_latex` API.

HTML output keeps selectable native MathML and embeds the same STIX Two Math font as
the documentation site, including in standalone offline notebooks. Math defaults to
21 px with padding and horizontal scrolling for long expressions. Set the CSS custom
property `--spenso-math-font-size` on the notebook or an output container to change
the base size; nested scripts retain their relative sizing.

Idenso transformations still return ordinary Symbolica expressions. Module-level
`spenso.formatted(expression)` or `spenso.as_tensor(expression)` provides tensor-aware display.
HTML, SVG, and rich-display functions accept `notation_source` as a trusted, complete
replacement for the bundled `notation.typ`, not a style fragment. Typst executes that source;
it must implement the expected notation interface and must not come from untrusted input.
Display customization stays outside Atom payloads; portable representation and math-label
declarations continue to travel with the expressions.

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

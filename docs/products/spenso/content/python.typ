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

`TensorExpression`, `Tensor`, and `TensorNetwork` share semantic display methods.
`DisplaySettings` controls the ports, Schoonschip, and call layouts, dimensions, parentheses,
commas, symbol scripts, and index/factor spacing. Positional calls such as `to_typst(True)`
and `formatted(True)` still request dimensions.

// docs-example: compile
```python
from symbolica.community.spenso import DisplaySettings, TensorExpression

trace = TensorExpression.gamma5(4).trace()
source = trace.to_typst(
    settings=DisplaySettings(show_dimensions=True, parentheses=False)
)
```

`to_typst` and `format_tensor` emit source using the ports layout. Schoonschip, call, and
custom-spacing settings require Tydenso's Typst notation layer; use HTML, SVG, or notebook
display for those settings. Source-only methods reject unsupported settings rather than
silently ignoring them.

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
LaTeX or text. `TensorNetwork.__str__` remains Graphviz DOT; `to_dot()` makes that intention
explicit. These display methods do not replace Symbolica's inherited `to_latex` API.

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

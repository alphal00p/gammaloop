#import "../../shared.typ": boundary, source-link

#let api = [
= Rust and Python APIs

== Rust package

The `idenso` crate exposes representation types and syntax macros together with several
rewrite families:

- `IndexTooling` covers canonicalization, conjugation, index wrapping, and dangling-index
  inspection for Symbolica atoms;
- `Cookable`, `CookSettings`, and the cook filters control reversible or flattening encodings;
- `dirac`, `color`, `epsilon`, and shorthand modules implement algebra-specific rewrites;
- `representations::initialize` installs the standard representation and tensor symbols.

Representation helper macros such as `bis!`, `cof!`, and `coad!` construct the symbolic forms
expected by Spenso and Idenso. The #link("reference/rust/")[Rust orientation] leads to the
revision-specific Rustdoc for their accepted forms, return types, feature gates, and source
locations. Native features such as `bincode`, `reference-cases`, and `wasm` are available only
when the matching Cargo feature is enabled. Python bindings belong to the separate `spynso3`
adapter.

== Python community module

#boundary("Part of Symbolica community", [
  Import the shared tensor API from `symbolica.community.tensor`. The `spynso3` adapter exposes
  Spenso tensors and Idenso algebra through one community module; `idenso` itself is a native
  Rust package, not a standalone PyPI distribution.
])

Install the published assembly with `python -m pip install --upgrade symbolica`, then verify it
with `python -c "import symbolica.community.tensor"`. There is no `pip install idenso` fallback.
For a source build, add `spynso3` to the external
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly and
register `spynso3::SpensoModule` through `SymbolicaCommunityModule` in that shared extension.
Building the native Idenso crate alone does not add the module to an installed Symbolica package.

The generated #link("reference/python/")[Python API] records exact signatures and defaults. Its
operations group naturally into four phases:

- construction: `Representation`, `TensorName`, and `TensorExpression` preserve typed slots;
- index preparation: `TensorExpression.list_dangling` and `wrap_indices` expose or isolate
  index interfaces;
- contraction: `TensorExpression.contract` contracts compatible tensor structures;
- algebra: `TensorExpression.simplify_algebra` selects gamma, color, epsilon, and contraction
  passes explicitly.

```python
from symbolica.community.tensor import Representation, TensorExpression, TensorName

minkowski = Representation.mink(4)
mu = minkowski("mu")
nu = minkowski("nu")
metric = TensorExpression.g(minkowski)
momentum = TensorName.vector("p")
expression = metric(mu, nu) * momentum(mu)

external_indices = expression.list_dangling()
reduced = expression.contract()
assert len(external_indices) == 1
assert len(reduced.list_dangling()) == 1
```

Idenso does not define a second parser syntax: the example constructs a Spenso-compatible
Symbolica expression and then applies one Idenso transformation. Keep transformations separate
while developing a pipeline so expression growth and convention changes remain observable.

For implementation details, start with
#source-link("crates/idenso/src/lib.rs", label: "the Rust API") and
#source-link("crates/spynso3/src/expression.rs", label: "the shared Python binding").
]

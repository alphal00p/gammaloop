#import "../../shared.typ": boundary, source-link

#let api = [
= Rust and Python APIs

== Rust package

The `idenso` crate exposes representation types and syntax macros together with several
rewrite families:

- `IndexTooling` covers canonicalization, conjugation, index wrapping, and dangling-index
  inspection for Symbolica atoms;
- `Cookable`, `CookSettings`, and the cook filters control reversible or flattening encodings;
- `SelectiveExpand` expands terms by recognized representation families;
- `dirac`, `color`, `epsilon`, and shorthand modules implement algebra-specific rewrites;
- `representations::initialize` installs the standard representation and tensor symbols.

Representation helper macros such as `bis!`, `cof!`, and `coad!` construct the symbolic forms
expected by Spenso and Idenso. The #link("reference/rust/")[Rust orientation] leads to the
revision-specific Rustdoc for their accepted forms, return types, feature gates, and source
locations. APIs behind `bincode`, `reference-cases`, `python`, and `python_stubgen` are available
only when the matching Cargo feature is enabled.

== Python community module

#boundary("Part of Symbolica community", [
  Import Idenso from `symbolica.community.idenso`. The `idenso` Cargo package supplies this
  community module when built with its Python feature; it is not a standalone PyPI distribution
  and is not included in the GammaLoop wheel merely because the crates share a repository.
])

Install the published assembly with `python -m pip install --upgrade symbolica`, then verify it
with `python -c "import symbolica.community.idenso"`. There is no `pip install idenso` fallback.
For a source build, add this crate to the external
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly with the
`python` feature and call `IdensoModule`'s `SymbolicaCommunityModule` registration when that
extension is assembled. Building the Rust crate alone does not add the community module to an
already installed Symbolica package.

The generated #link("reference/python/")[Python API] records exact signatures and defaults. Its
operations group naturally into four phases:

- setup: importing the community module registers its symbols;
- expansion: `expand_bis`, `expand_mink`, `expand_mink_bis`, `expand_metrics`, and
  `expand_color`;
- index preparation: `wrap_indices`, `wrap_dummies`, `list_dangling`, `cook_indices`, and
  `cook_function`;
- algebra: `dirac_adjoint`, `simplify_gamma`, `to_dots`, `simplify_metrics`, and
  `simplify_color`.

```python
from symbolica.community.idenso import list_dangling, simplify_metrics
from symbolica.community.spenso import Representation, TensorName

minkowski = Representation.mink(4)
mu = minkowski("mu")
nu = minkowski("nu")
metric = TensorName.g()
momentum = TensorName("p")
expression = metric(mu, nu) * momentum(mu)

external_indices = list_dangling(expression)
reduced = simplify_metrics(expression)
assert len(external_indices) == 1
assert len(list_dangling(reduced)) == 1
```

Idenso does not define a second parser syntax: the example constructs a Spenso-compatible
Symbolica expression and then applies one Idenso transformation. Keep transformations separate
while developing a pipeline so expression growth and convention changes remain observable.

For implementation details, start with
#source-link("crates/idenso/src/lib.rs", label: "the Rust API") and
#source-link("crates/idenso/src/python.rs", label: "the Python binding").
]

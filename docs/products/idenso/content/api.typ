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
locations. APIs behind `bincode` and `reference-cases` are available only when the matching Cargo
feature is enabled. Python bindings belong to `spynso3`.

== Python community module

#boundary("Part of Symbolica community", [
  Import `TensorExpression` and algebra settings from `symbolica.community.spenso`.
  The `spynso3` Cargo package supplies the unified Python bindings for Spenso and Idenso;
  Idenso remains the Rust algebra implementation. There is no separate `idenso` Python module.
])

Install a Symbolica community build with the unified `spynso3` bindings, then verify it with
`python -c "import symbolica.community.spenso"`. There is no `pip install idenso` fallback.
Source embedders add `spynso3` to their
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly and
register `SpensoModule` through `SymbolicaCommunityModule`. Building the Rust crate alone does
not add the community module to an already installed Symbolica package.

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
from symbolica.community.spenso import Representation, TensorExpression, TensorName

minkowski = Representation.mink(4)
mu = minkowski("mu")
nu = minkowski("nu")
metric = TensorExpression.g(minkowski)
momentum = TensorName.vector("p")
expression = metric(mu, nu) * momentum(mu)

external_indices = expression.list_dangling()
reduced = expression.schoonschip_net()
assert len(external_indices) == 1
assert len(reduced.list_dangling()) == 1
assert reduced == momentum(nu)
```

Use `schoonschip_net()` to execute typed `bracket` contractions. Focused simplifiers such as
`simplify_metrics()` rewrite ordinary indexed products.

Idenso does not define a second parser syntax: the example constructs a Spenso-compatible
`TensorExpression` and then applies one Idenso transformation through its methods. Keep transformations separate
while developing a pipeline so expression growth and convention changes remain observable.

For implementation details, start with
#source-link("crates/idenso/src/lib.rs", label: "the Rust API") and
#source-link("crates/spynso3/src/expression.rs", label: "the Python binding").
]

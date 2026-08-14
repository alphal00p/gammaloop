#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python API boundaries

#catalog-contract(
  rust-scope: "idenso",
  python-scope: "symbolica.community.idenso",
)

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
expected by Spenso and Idenso. Exact macro grammar and rewrite return types come from the
generated Rust catalog. The optional `bincode`, `reference-cases`, `python`, and
`python_stubgen` surfaces must be labeled by feature.

== Python community module

#boundary("Part of Symbolica community", [
  Import Idenso from `symbolica.community.idenso`. The `idenso` Cargo package supplies this
  community module when built with its Python feature; it is not a standalone PyPI distribution
  and is not included in the GammaLoop wheel merely because the crates share a repository.
])

Install a Symbolica Python distribution only if its release manifest says it includes Idenso,
then verify the actual assembly with
`python -c "import symbolica.community.idenso"`. There is no `pip install idenso` fallback.
For a source build, add this crate to the external
#link("https://github.com/benruijl/symbolica-community")[symbolica-community] assembly with the
`python` feature and call `IdensoModule`'s `SymbolicaCommunityModule` registration when that
extension is assembled. This documentation build validates the registration function but does
not modify an installed Symbolica package.

The Python API operates on Symbolica expressions and groups naturally into four phases:

- setup: `initialize`;
- expansion: `expand_bis`, `expand_mink`, `expand_mink_bis`, `expand_metrics`, and
  `expand_color`;
- index preparation: `wrap_indices`, `wrap_dummies`, `list_dangling`, `cook_indices`, and
  `cook_function`;
- algebra: `dirac_adjoint`, `simplify_gamma`, `to_dots`, `simplify_metrics`, and
  `simplify_color`.

```python
from symbolica.community.idenso import (
    initialize,
    list_dangling,
    simplify_metrics,
)

initialize()
# `expr` is a Symbolica Expression written with Spenso-compatible tensor forms.
external_indices = list_dangling(expr)
reduced = simplify_metrics(expr)
```

The example deliberately leaves expression construction to Symbolica and Spenso: Idenso does
not define a second parser syntax. Apply one transformation at a time while developing a
pipeline so expression growth and convention changes remain observable.

Python docstrings, signatures, and module names come from the binding descriptor and are
rendered as plain text. Source starting points are
#source-link("crates/idenso/src/lib.rs", label: "the Rust API") and
#source-link("crates/idenso/src/python.rs", label: "the Python binding").
]

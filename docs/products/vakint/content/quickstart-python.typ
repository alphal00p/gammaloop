#import "../../shared.typ": callout

#let quickstart-python = [
= Using Vakint from Python

Vakint is available as `symbolica.community.hepkit.vakint` in a Symbolica community
assembly built with this checkout. This matching-only workflow canonicalizes one
loop with arbitrary input labels and invokes no external evaluation tool.

== Verify the installed module

Install the combined host wheel built with this checkout, then verify its public path:

// docs-example: syntax
```sh
python -c "import symbolica.community.hepkit.vakint as vakint; print(vakint.__name__)"
```

For source assemblies, `VakintWrapper::get_name()` returns `hepkit_vakint`.
Register its native module as `symbolica.community.hepkit_vakint_native` and copy
`crates/vakint/python/symbolica/community/hepkit/vakint/` into the host's matching
Python package directory. The host's `symbolica.community.hepkit` must be a package;
the wrapper imports the native exports and calls `initialize_module()`.

There is no separate `vakint` Python wheel. Follow
#link("https://symbolica.io/docs/get_started.html")[Symbolica's installation and license terms].

== Canonicalize a one-loop integral

Save this as `vakint_quickstart.py`:

// docs-example: compile vakint-community-quickstart
```python
from symbolica import E
from symbolica.community.hepkit.vakint import Vakint

engine = Vakint(evaluation_order=[])
integral = E(
    "topo(prop(18,edge(7,7),k(99),muvsq,1))",
    default_namespace="vakint",
)
canonical = engine.to_canonical(integral, short_form=True)

assert "I1L" in str(canonical)
print(canonical)
```

Run `python vakint_quickstart.py`. The arbitrary propagator, edge, and momentum labels are
normalized to Vakint's one-loop topology. An empty `evaluation_order` is intentional: it prevents
this first use from probing FORM, MATAD, FMFT, or pySecDec.

#callout("Matching and evaluation are separate choices", [
  Canonicalization answers which supported topology an expression represents. Tensor reduction
  and numerical or analytic evaluation add backend requirements and normalization choices.
  Enable them only after the canonical form is understood.
])

Use the #link("quickstart/rust/")[Rust guide] for the native matching API, or continue with
the #link("guides/evaluation/")[matching and evaluation guide] before enabling a backend.
]

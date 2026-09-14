#import "../../shared.typ": callout

#let quickstart-python = [
= Using Vakint from Python

Symbolica `2.2.0` bundles Vakint's community module. This matching-only workflow canonicalizes one
loop with arbitrary input labels and invokes no external evaluation tool.

== Install and verify the module

// docs-example: syntax
```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "symbolica==2.2.0"
python -c "import symbolica.community.vakint as vakint; print(vakint.__name__)"
```

There is no separate `vakint` Python wheel. Follow
#link("https://symbolica.io/docs/get_started.html")[Symbolica's installation and license terms].

== Canonicalize a one-loop integral

Save this as `vakint_quickstart.py`:

// docs-example: compile vakint-community-quickstart
```python
from symbolica import E
from symbolica.community.vakint import Vakint

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

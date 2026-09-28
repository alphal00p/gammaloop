#import "../../shared.typ": callout

#let quickstart-python = [
= Using Idenso from Python

The `spynso3` bindings expose Spenso and Idenso through one Symbolica community module. This example contracts one
four-dimensional Minkowski metric with a vector and verifies the exact remaining expression.

== Install and verify the modules

// docs-example: syntax
```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install symbolica
python -c "import symbolica.community.spenso"
```

There is no standalone `idenso` Python wheel. Import the community module before constructing
expressions, and follow
#link("https://symbolica.io/docs/get_started.html")[Symbolica's installation and license terms].

== Simplify one metric contraction

Save this as `idenso_quickstart.py`:

// docs-example: compile idenso-community-quickstart
```python
from symbolica.community.spenso import Representation, TensorExpression, TensorName

rep = Representation.mink(4)
mu, nu = rep("mu"), rep("nu")
g = TensorExpression.g(rep)
q = TensorName.vector("q")

expression = g(mu, nu) * q(mu)
reduced = expression.contract().to_expression()

assert reduced == q(nu)
assert len(reduced.list_dangling()) == 1
print(reduced)
```

Typed composition represents the contraction with a `bracket` shorthand. `contract()`
uses the shared symbolic contractor for bracketed and ordinary indexed products. It returns
an `AliasedTensorExpression`; `to_expression()` resolves its literal definitions without
distributing scalar coefficients. The #link("guides/algebra/")[algebra guide] explains dummy scopes.

Run `python idenso_quickstart.py`. Success means the explicit metric disappears and `q(nu)`
remains. The structural equality check is stronger than comparing printed text, whose formatting
may change between Symbolica releases.

#callout("Make rewrite stages observable", [
  Check dangling indices before and after each pass. When multiplying expressions built with
  independent dummy-index namespaces, wrap those namespaces before simplification so equal
  printed names do not create an accidental contraction.
])

Continue with the #link("tutorial/")[controlled identity tutorial] for a verification workflow,
then use the #link("guides/algebra/")[algebra guide] to compose metric, Dirac, color, and cooking
stages deliberately.
]

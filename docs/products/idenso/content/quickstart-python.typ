#import "../../shared.typ": callout

#let quickstart-python = [
= Using Idenso from Python

Symbolica `2.2.0` bundles Idenso and Spenso as native community modules. This example contracts one
four-dimensional Minkowski metric with a vector and verifies the exact remaining expression.

== Install and verify the modules

// docs-example: syntax
```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "symbolica==2.2.0"
python -c "import symbolica.community.idenso; import symbolica.community.spenso"
```

There is no standalone `idenso` Python wheel. Import the community module before constructing
expressions, and follow
#link("https://symbolica.io/docs/get_started.html")[Symbolica's installation and license terms].

== Simplify one metric contraction

Save this as `idenso_quickstart.py`:

// docs-example: compile idenso-community-quickstart
```python
from symbolica.community.idenso import initialize, list_dangling, simplify_metrics
from symbolica.community.spenso import Representation, TensorName

initialize()
rep = Representation.mink(4)
mu, nu = rep("mu"), rep("nu")
g, q = TensorName.g(), TensorName("q")

expression = g(mu, nu) * q(mu)
reduced = simplify_metrics(expression)

assert reduced == q(nu).to_expression()
assert len(list_dangling(reduced)) == 1
print(reduced)
```

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

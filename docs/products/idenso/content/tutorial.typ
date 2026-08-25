#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

This tutorial uses Idenso's Python community module to contract one four-dimensional Minkowski
metric tensor with a vector. It is intentionally small: the point is to establish the registered
tensor syntax and observe one algebra pass before composing a larger Dirac or color pipeline.

== Prerequisites

Idenso is mounted at `symbolica.community.idenso`; it is not installed by `pip install idenso`.
Install the published Symbolica distribution, which bundles the Idenso and Spenso community
modules, and verify the actual environment first:

// docs-example: syntax
```sh
python -m pip install --upgrade symbolica
python -c "import symbolica.community.idenso; import symbolica.community.spenso"
```

If that fails, check the installed Symbolica version before continuing. Generating a `.pyi`
file does not mount the native module. Source embedders can build a custom
#link("https://github.com/symbolica-dev/symbolica-community")[community assembly].

== Simplify one metric contraction

Save the following as `metric_first.py`:

// docs-example: compile idenso-controlled-identity
```python
from symbolica.community.idenso import list_dangling, simplify_metrics
from symbolica.community.spenso import Representation, TensorName

rep = Representation.mink(4)
mu = rep("mu")
nu = rep("nu")

g = TensorName.g()
q = TensorName("q")
expression = g(mu, nu) * q(mu)

free_before = list_dangling(expression)
reduced = simplify_metrics(expression)
free_after = list_dangling(reduced)

assert len(free_before) == 1
assert len(free_after) == 1
print("reduced:", reduced)
```

Run `python metric_first.py`. Success means both rank assertions pass, the metric is removed from
the reduced expression, and the result is a rank-one expression carrying `nu`. The rewrite is
abstract: the registered Minkowski metric supplies index compatibility, but this example does
not substitute a numerical signature. Printed output can vary with the installed Symbolica
version; inspect the expression structure instead of comparing exact text.

#callout("Verification scope and cost", [
  The docs harness compiles this Python source without importing native modules; it syntax-checks
  the environment probe without running it. Execute both steps in a provisioned community-module
  environment to verify the rank-one invariant. The script should take seconds after startup;
  building or installing Symbolica community modules is the expensive prerequisite and may take
  substantially longer.
])

#callout("Import before constructing expressions", [
  Importing `symbolica.community.idenso` registers Idenso's representation and tensor symbols.
  Construct the Spenso-compatible expression after that import so Idenso's matchers see the
  intended slots and tensor names.
])

== Grow this into an algebra pipeline

Keep transformations observable while developing:

- call `list_dangling` before and after a pass to catch accidental contractions;
- use `wrap_dummies` before multiplying expressions built in independent index namespaces;
- expand only the Minkowski, bispinor, or color sector needed by the next pass;
- apply `simplify_metrics`, `simplify_gamma`, and `simplify_color` as distinct phases;
- canonicalize only after the physics-specific identities required by the calculation are
  explicit.

== Troubleshooting and next steps

- `ModuleNotFoundError` for `symbolica.community.idenso` means the installed Symbolica build
  does not include the native module; a `.pyi` type stub or source checkout alone is not enough.
- If the metric remains unchanged, construct its slots with the same `Representation` and
  abstract index objects as the vector. Plain Symbolica functions do not automatically carry
  Spenso tensor structure.
- If a free index disappears unexpectedly, inspect repeated names and use `wrap_dummies` before
  combining separately constructed expressions.
- Continue with the syntax-and-algebra manual, then use the Python and Rust API references for
  signatures and feature gates.
]

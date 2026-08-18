#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

This tutorial uses Idenso's Python community module to contract one metric tensor with a
vector. It is intentionally small: the point is to establish the registered tensor syntax and
observe one algebra pass before composing a larger Dirac or color pipeline.

== Prerequisites

Idenso is mounted at `symbolica.community.idenso`; it is not installed by `pip install idenso`.
Use a Symbolica Python distribution whose release includes the Idenso and Spenso community
modules, or build those modules through the `symbolica-community` assembly. Verify the actual
environment first:

The assembly and installation instructions live in
#link("https://github.com/benruijl/symbolica-community")[`symbolica-community`].

// docs-example: syntax
```sh
python -c "import symbolica.community.idenso; import symbolica.community.spenso"
```

If that fails, install or rebuild the community assembly before continuing. Generating a `.pyi`
file does not mount the native module.

== Simplify one metric contraction

Save the following as `metric_first.py`:

// docs-example: compile idenso-controlled-identity
```python
from symbolica.community.idenso import initialize, list_dangling, simplify_metrics
from symbolica.community.spenso import Representation, TensorName

initialize()
rep = Representation.euc(3)
mu = rep("mu")
nu = rep("nu")

g = TensorName.g()
q = TensorName("q")
expression = g(mu, nu) * q(mu)

print("free before:", list_dangling(expression))
reduced = simplify_metrics(expression)
print("reduced:", reduced)
print("free after:", list_dangling(reduced))
```

Run `python metric_first.py`. Success means the metric is removed from the reduced expression,
the contracted `mu` no longer appears as a free index, and the result is a rank-one expression
carrying `nu`. Printed output can vary with the installed Symbolica version; inspect the
expression structure instead of comparing exact text.

#callout("Verification scope and cost", [
  The docs harness compiles this Python source without importing native modules; it syntax-checks
  the environment probe without running it. Execute both steps in a provisioned community-module
  environment to verify the rank-one invariant. The script should take seconds after startup;
  building or installing Symbolica community modules is the expensive prerequisite and may take
  substantially longer.
])

#callout("Initialize before parsing or rewriting", [
  `initialize()` installs Idenso's representation and tensor symbols. Calling it explicitly at
  the start of a standalone script makes registration order clear and keeps expressions in the
  Spenso-compatible form that Idenso's matchers expect.
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

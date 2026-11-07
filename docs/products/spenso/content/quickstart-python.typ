#import "../../shared.typ": callout

#let quickstart-python = [
= Using Spenso from Python

Spenso is available as `symbolica.community.tensor` in a Symbolica Community assembly.
This example uses the checkout's combined Symbolica `3.0.1` host, which includes Spenso,
Linnet and FeynKit, to construct a typed two-dimensional tensor.

== Install and verify the module

From the repository root, create an isolated environment and install the locked Community
host with the native development command:

// docs-example: syntax
```sh
python -m venv .venv
. .venv/bin/activate
nix develop -c just python-develop community --python .venv/bin/python
python -c "import symbolica.community.tensor as spenso; print(spenso.__name__)"
```

There is no separate `spenso` Python wheel. The installed `symbolica` distribution owns the
native `symbolica.community.tensor` module; the version number alone does not establish which
Community extensions an external wheel includes. Keep the generated `.venv/.nix-python/symbolica`
root while using this Nix-built host: its native libraries remain in the Nix store. Symbolica's
#link("https://symbolica.io/docs/get_started.html")[installation and license terms] apply; its
free restricted mode is enough for this single-process example.

== Construct a dense tensor

Save this as `spenso_quickstart.py`:

// docs-example: compile spenso-community-quickstart
```python
from symbolica.community.tensor import Representation, Tensor, TensorName

rep = Representation.euc(2)
matrix = Tensor.dense(
    TensorName("quickstart::matrix")(rep("i"), rep("j")),
    [1.0, 0.0, 0.0, 1.0],
)

assert matrix[0, 0] == 1.0
assert matrix[1, 1] == 1.0
matrix.to_sparse()
print(matrix)
```

Run `python spenso_quickstart.py`. `TensorName` gives the data a named, representation-aware
structure, while the four row-major values provide its concrete storage. Converting to sparse
storage changes only that storage; it does not change the slots or re-index the tensor.

#callout("Structure comes before storage", [
  Equal dimensions alone do not make two tensor slots compatible. Their representations,
  variance, and abstract indices remain part of the value. Diagnose those properties before
  changing dense or sparse storage.
])

Continue with the #link("guides/python/")[Python tensor-workflow guide], or use the
#link("quickstart/rust/")[Rust guide] to perform a checked native contraction.
]

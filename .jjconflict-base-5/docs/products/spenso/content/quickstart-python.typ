#import "../../shared.typ": callout

#let quickstart-python = [
= Using Spenso from Python

The shortest path into Spenso is the community module bundled with Symbolica `2.2.0`. This
example creates a typed two-dimensional tensor and checks its diagonal entries; it needs no
GammaLoop checkout.

== Install and verify the module

Create an isolated environment, install the version used by this manual, and verify the community
import:

// docs-example: syntax
```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "symbolica==2.2.0"
python -c "import symbolica.community.spenso as spenso; print(spenso.__name__)"
```

There is no separate `spenso` Python wheel. The installed `symbolica` distribution owns the
native `symbolica.community.spenso` module. Symbolica's
#link("https://symbolica.io/docs/get_started.html")[installation and license terms] apply; its
free restricted mode is enough for this single-process example.

== Construct a dense tensor

Save this as `spenso_quickstart.py`:

// docs-example: compile spenso-community-quickstart
```python
from symbolica.community.spenso import Representation, Tensor, TensorIndices

rep = Representation.euc(2)
matrix = Tensor.dense(
    TensorIndices(rep("i"), rep("j")),
    [1.0, 0.0, 0.0, 1.0],
)

assert matrix[0, 0] == 1.0
assert matrix[1, 1] == 1.0
matrix.to_sparse()
print(matrix)
```

Run `python spenso_quickstart.py`. `TensorIndices` gives the data a representation-aware
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

"""Exercise constructor selectors against an installed community extension."""

import importlib
import sys

from symbolica import E, S
from symbolica.community.tensor import Representation, TensorName

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
D, mu = S("selector_test::D", "selector_test::mu")
mink, dot = S("spenso::mink", "spenso::dot")
lorentz = Representation.mink(D)
name = TensorName.vector("selector_test::k")
other = TensorName.vector("other_namespace::k")
k, p = name.to_expression(), other.to_expression()
kc, pc = k(0, mink(D)), p(0, mink(D))
ki, pi = k(0, mink(D, mu)), p(0, mink(D, mu))

# TensorName and bare Expression select the entire head, preserving namespaces.
for head in (name, k):
    reducer = fk.TensorReducer(D, integrated=[head])
    assert reducer.reduce(ki) == E("0")
    assert reducer.reduce(k(1, mink(D, mu))) == E("0")
    assert reducer.reduce(pi) == pi

# Expression and TensorExpression select an exact compact vector, including args.
for vector in (kc, name(0, lorentz)):
    reducer = fk.TensorReducer(D, integrated=(vector,))
    assert reducer.reduce(ki) == E("0")
    assert reducer.reduce(k(1, mink(D, mu))) == k(1, mink(D, mu))
    assert reducer.reduce(pi) == pi

# Both external representations produce the same longitudinal projection.
expected = dot(kc, pc) / dot(pc, pc) * pi
for external in (pc, other(0, lorentz)):
    reducer = fk.TensorReducer(D, integrated=[name(0, lorentz)], external=[external])
    assert (reducer.reduce(ki) - expected).together() == E("0")
    assert reducer.reduce(pi) == pi

# An external basis entry takes precedence over selecting its whole head.
reducer = fk.TensorReducer(D, integrated=[name, other], external=[other(0, lorentz)])
assert (reducer.reduce(ki) - expected).together() == E("0")
assert reducer.reduce(pi) == pi

for invalid in (E("0"), k + 1, k**2):
    try:
        fk.TensorReducer(D, integrated=[invalid])
    except ValueError:
        pass
    else:
        raise AssertionError(f"accepted invalid selector: {invalid}")
for keyword, invalid in (("integrated", "k"), ("external", "p"), ("external", other)):
    try:
        fk.TensorReducer(D, **{keyword: [invalid]})
    except TypeError:
        pass
    else:
        raise AssertionError(f"accepted invalid {keyword} entry: {invalid}")

empty = fk.TensorReducer(D)
try:
    empty.reduce(ki)
except fk.TensorReductionError:
    pass
else:
    raise AssertionError("reduced without selecting integrated momenta")
for removed in (
    "with_integrated_head",
    "with_integrated_vector",
    "with_external_vector",
):
    assert not hasattr(empty, removed)
print("Tensor constructor selectors passed")

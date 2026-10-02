"""Keep four-loop gluon tensor projection factored in an installed host."""

import importlib
import sys
from time import perf_counter

from symbolica import E
from symbolica.community.tensor import TensorExpression

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
diagrams = (
    fk.Model.qcd()
    .process(["g"], ["g"], particle_veto=["u", "c", "s", "t", "b"])
    .generate_diagrams(loops=4)
    .diagrams
)
gluons = [
    diagram
    for diagram in diagrams
    if len(diagram.internal_edges) == 11
    and all(edge.particle_name == "g" for edge in diagram.internal_edges)
]
assert gluons
for diagram in (gluons[0], gluons[len(gluons) // 2], gluons[-1]):
    original = diagram.to_json()
    projected = diagram.numerator_expression() * diagram.projector_expression()
    start = perf_counter()
    reduced = diagram.tensor_reduce(E("D"))
    elapsed = perf_counter() - start
    assert isinstance(reduced, TensorExpression)
    assert reduced.structure == projected.structure
    assert diagram.to_json() == original
    # The regression used to distribute eight six-term vertices before
    # projection. Inspect the retained expression without expanding it.
    assert len(str(reduced.to_expression())) < 200_000
    print(f"{diagram.name}: {elapsed:.3f}s")

"""Check diagram UV expansion in an installed FeynKit or HEP host."""

import importlib
import sys
from pathlib import Path

from symbolica import E, S
from symbolica.community.spenso import TensorExpression

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
bubble = (
    model.process(["scalar_0"], ["scalar_0"], vertex_allow=["V_3_SCALAR_000"])
    .generate_diagrams(loops=1, max_vertices=2, allow_self_loops=False)
    .diagrams[0]
)
mass = S("uv_test::mUV", is_scalar=True)
original = bubble.to_json()
expanded = bubble.uv_expansion(mass)
assert isinstance(expanded, TensorExpression)
assert expanded.is_scalar
assert expanded != 0
assert bubble.uv_counterterm(mass) == -expanded
assert bubble.to_json() == original
region = bubble.filter(
    edge=lambda e: e.data.id in {x.id for x in bubble.internal_edges}
)
assert region.uv_expansion(mass) == expanded
assert region.uv_counterterm(mass) == -expanded
assert bubble.uv_expansion(mass, dimension=2) == 0
assert bubble.uv_expansion(mass, numerator=0) == 0
assert bubble.subgraph().uv_expansion(mass) == 0
tree = bubble.filter(edge=lambda e: e.data.id == bubble.internal_edges[0].id)
assert tree.uv_expansion(mass) == 0

# Signed powers share denominator_expression's edge IDs and selection semantics.
internal_ids = [edge.id for edge in bubble.internal_edges]
for powers in ({internal_ids[0]: 2}, {internal_ids[0]: 0}, {internal_ids[0]: -1}):
    powered = bubble.uv_expansion(mass, dimension=6, numerator=1, edge_powers=powers)
    assert (
        region.uv_expansion(mass, dimension=6, numerator=1, edge_powers=powers)
        == powered
    )
    assert (
        bubble.uv_counterterm(mass, dimension=6, numerator=1, edge_powers=powers)
        == -powered
    )
    assert (
        region.uv_counterterm(mass, dimension=6, numerator=1, edge_powers=powers)
        == -powered
    )
assert bubble.uv_expansion(mass, numerator=1, edge_powers={internal_ids[0]: 2}) == 0
assert bubble.uv_expansion(mass, edge_powers={10**6: 2}) == expanded
for operation in (
    bubble.denominator_expression,
    bubble.uv_expansion,
    bubble.uv_counterterm,
):
    arguments = () if operation.__name__ == "denominator_expression" else (mass,)
    try:
        operation(*arguments, edge_powers={-1: 2})
    except OverflowError:
        pass
    else:
        raise AssertionError("negative edge ID accepted")

# An external-only numerator must stay soft even when its polynomial degree is high.
basis = region.momentum_basis()
external = next(e for e in basis.external_edges if e not in basis.dependent_externals)
soft = E(f"gammalooprs::Q({external},spenso::mink(4,uv_test::mu))")
tensor = bubble.uv_expansion(mass, numerator=soft * bubble.numerator_expression())
assert isinstance(tensor, TensorExpression)
assert tensor.rank == 1
assert tensor == soft * expanded
assert bubble.uv_expansion(mass, numerator=soft**4, dimension=2) == 0

# Filters carry diagram ownership, even for an identical serialized topology.
copy = fk.FeynmanDiagram.from_json(model, original)
try:
    copy.subgraph(region)
except (ValueError, TypeError):
    pass
else:
    raise AssertionError("foreign filter accepted")
try:
    bubble.uv_expansion(mass, dimension=0)
except fk.DiagramError:
    pass
else:
    raise AssertionError("nonpositive dimension accepted")

triangle = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_diagrams(loops=1, max_vertices=3, allow_self_loops=False)
    .diagrams
)
triangle = next(
    d
    for d in triangle
    if len({e.source if e.source is not None else e.target for e in d.external_edges})
    == 3
)
assert triangle.uv_counterterm(mass) == 0
assert triangle.uv_counterterm(mass, dimension=6) != 0

# The HEP example's gluon, ghost and massive-quark bubbles retain their open
# Lorentz/color interface; callers can contract a projector after expansion.
sm = fk.Model(Path(__file__).parents[3] / "assets/models/json/sm/sm.json")
diagrams = (
    sm.process(["g"], ["g"], particle_veto=["c", "t", "s", "u", "d"])
    .generate_diagrams(
        loops=1,
        max_vertices=2,
        coupling_orders={"QCD": 2, "QED": 0},
        allow_self_loops=False,
    )
    .diagrams
)
assert len(diagrams) >= 3
for diagram in diagrams:
    counterterm = diagram.uv_counterterm(mass)
    assert isinstance(counterterm, TensorExpression)
    assert counterterm != 0
    assert counterterm.rank == 4
    assert set(counterterm.list_dangling()) == set(
        diagram.numerator_expression().list_dangling()
    )

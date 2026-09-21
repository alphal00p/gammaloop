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
bubble = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0"],
    loops=1,
    max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=False,
).diagrams[0]
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
assert bubble.uv_expansion(mass, subgraph=region) == expanded
assert bubble.uv_counterterm(mass, subgraph=region) == -expanded
assert bubble.uv_expansion(mass, dimension=2) == 0
assert bubble.uv_expansion(mass, numerator=0) == 0
assert bubble.uv_expansion(mass, subgraph=bubble.to_linnet().empty_subgraph()) == 0
tree = bubble.filter(edge=lambda e: e.data.id == bubble.internal_edges[0].id)
assert bubble.uv_expansion(mass, subgraph=tree) == 0

# An external-only numerator must stay soft even when its polynomial degree is high.
basis = bubble.momentum_basis(subgraph=region)
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
    copy.uv_counterterm(mass, subgraph=region)
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

triangle = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=1,
    max_vertices=3,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=False,
).diagrams
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
diagrams = sm.generate_diagrams(
    ["g"],
    ["g"],
    loops=1,
    max_vertices=2,
    coupling_orders={"QCD": 2, "QED": 0},
    particle_veto=["c", "t", "s", "u", "d"],
    allow_self_loops=False,
).diagrams
assert len(diagrams) >= 3
for diagram in diagrams:
    counterterm = diagram.uv_counterterm(mass)
    assert isinstance(counterterm, TensorExpression)
    assert counterterm != 0
    assert counterterm.rank == 4
    assert set(counterterm.list_dangling()) == set(
        diagram.numerator_expression().list_dangling()
    )

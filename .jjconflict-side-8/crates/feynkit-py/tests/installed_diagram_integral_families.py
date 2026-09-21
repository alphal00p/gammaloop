"""Generated diagrams use shared propagators, routing, and family algorithms."""

from pathlib import Path

from symbolica import E, S
from symbolica.community import feynkit as fk

model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0"],
    loops=1,
    max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=False,
).diagrams
assert len(diagrams) == 1
s, x, y = S("s", "x", "y")
diagram = diagrams[0]
family = diagram.integral_family()
assert len(family.loop_momenta) == 1
assert len(family.external_momenta) == 1
assert len(family.denominators) == len(diagram.internal_edges) == 2
p = family.external_momenta[0]
kin = family.kinematics.with_scalar_product(p, p, s)
reference = diagram.integral_family(kinematics=kin)
for basis in diagram.loop_momentum_bases():
    routed = diagram.with_loop_momentum_edges(basis.loop_edges)
    family = routed.integral_family(kinematics=kin)
    assert family.is_complete and family.is_independent
    U, F = family.symanzik([x, y])
    assert U == x + y
    assert (F + s * x * y).expand() == E("0")
    assert family.find_mapping(reference) is not None
    assert family.scaleless_scaling([x, y]) is None

onshell = diagram.integral_family(kinematics=kin.with_scalar_product(p, p, E("0")))
assert onshell.scaleless_scaling([x, y]) is not None
four_dimensional = diagram.integral_family(kinematics=fk.Kinematics())
assert four_dimensional.kinematics.dimension == E("4")

contact = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=0,
    max_vertices=1,
    vertex_allow=["V_3_SCALAR_000"],
).diagrams[0]
try:
    contact.integral_family()
except fk.DiagramError:
    pass
else:
    raise AssertionError("tree diagram was accepted as a loop-integral family")
print("Generated diagram families, routing invariance and scaleless limits passed")

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

sunrise = fk.FeynmanDiagram.from_dot(
    fk.Model.phi4(),
    """digraph sunrise {
        ext [style=invis];
        ext -> a [particle="phi"];
        a -> b [particle="phi", lmb_id=0];
        a -> b [particle="phi", lmb_id=1];
        a -> b [particle="phi"];
        b -> ext [particle="phi"];
    }""",
)
raw = sunrise.propagator_family()
automatic = sunrise.integral_family()
assert len(raw.denominators) == 3
assert not raw.is_complete
assert automatic.is_complete and automatic.is_independent
assert automatic.denominators[:3] == raw.denominators
assert len(automatic.denominators) == 5
products = [
    raw.kinematics.scalar_product(k, raw.external_momenta[0]) for k in raw.loop_momenta
]
preferred = sunrise.integral_family(products)
assert (
    preferred.denominators
    == fk.IntegralFamily.from_diagram(sunrise, products).denominators
)
assert preferred.denominators == raw.denominators + products
partial = sunrise.integral_family(
    independent_dot_products=[raw.denominators[0], products[1]]
)
assert partial.denominators[3] == products[1]
assert partial.is_complete and partial.is_independent
assert fk.IntegralFamily.from_diagram(
    sunrise, kinematics=fk.Kinematics()
).kinematics.dimension == E("4")
try:
    fk.IntegralFamily.from_diagram(sunrise, products + [products[0] ** 2])
except fk.DiagramError:
    pass
else:
    raise AssertionError("nonlinear auxiliary product was accepted")

contact = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=0,
    max_vertices=1,
    vertex_allow=["V_3_SCALAR_000"],
).diagrams[0]
for construct in (
    lambda: contact.integral_family(),
    lambda: contact.propagator_family(),
    lambda: fk.IntegralFamily.from_diagram(contact),
):
    try:
        construct()
    except fk.DiagramError:
        pass
    else:
        raise AssertionError("tree diagram was accepted as a loop-integral family")
print("Generated diagram families, routing invariance and scaleless limits passed")

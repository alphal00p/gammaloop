"""Generated QED squared matrix element, including graph signs and propagators.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-MuAmu
Polarized production remains a separate benchmark.
"""

from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community import hep as fk
from symbolica.community.spenso import TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
# A sewn tree x tree graph has one loop and its two muon lines form the cut.
process = model.process([11, -11], [13, -13])
result = process.with_filters(vertex_allow=["V_98", "V_99"]).generate_cross_section(
    loops=1,
    max_vertices=4,
    maximum_bridges=None,
    numerator_grouping=None,
    progress=None,
)
assert len(result.diagrams) == 1
diagram = result.diagrams[0]
assert len(diagram.cuts) == 1
cut = diagram.cuts[0]
assert sorted(p.pdg_code for p in cut.particles) == [-13, 13]
assert all(side.loop_count == 0 for side in (cut.left, cut.right))

Q, P, K = S("gammalooprs::Q", "gammalooprs::P", "gammalooprs::K")
p1, p2, k1, k2 = P(0), P(1), S("k1"), S("k2")
index = S("index_")
s, t, u, me, mm, e = S("s", "t", "u", "UFO::Me", "UFO::MM", "UFO::ee")
kin = fk.Kinematics.mandelstam(
    [p1, p2, k1, k2], [me**2, me**2, mm**2, mm**2], [s, t, u]
)
# Select the mu- line explicitly; its cut orientation fixes the physical momentum.
muon_edge = next(edge for edge, p in zip(cut.edges, cut.particles) if p.pdg_code == 13)
diagram = diagram.with_loop_momentum_edges([muon_edge.id])
cut = diagram.cuts[0]
assert diagram.loop_momentum_basis.loop_edges == [muon_edge.id]
orientation = cut.orientations[muon_edge.id]

projector = diagram.projector_expression()
for edge in diagram.external_edges:
    particle = model.particle_by_pdg(edge.particle_pdg)
    projector = particle.sum_spins(projector, Q(edge.id), edge=edge.id, average=True)
    assert (
        particle.sum_spins(projector, Q(edge.id), edge=edge.id, average=True)
        == projector
    )

numerator = model.expand_couplings(
    diagram.numerator_expression().to_expression() * projector
)
contracted = (
    TensorExpression(numerator)
    .simplify_gamma()
    .to_expression()
    .to_dots()
    .to_expression()
)
contracted = kin.apply(
    diagram.loop_momentum_basis.route_expression(contracted).replace(
        K(0, index), orientation * k1(index)
    )
)

# Cut propagators belong to the phase-space measure, not the squared amplitude.
denominator = diagram.denominator_expression(
    edge_powers={edge.id: 0 for edge in cut.edges}, dimension=4, in_lmb=True
).to_expression()
denom = S("gammalooprs::denom")
a, b, c, inverse = S("a_", "b_", "c_", "inverse_")
denominator = kin.apply(
    denominator.replace(denom(a, b, c, inverse), inverse).replace(
        K(0, index), orientation * k1(index)
    )
)
assert (denominator - s**2).expand() == E("0")
factor = diagram.overall_factor_expression(evaluate=True)
assert factor == E("-1")
squared = factor * diagram.numerator_prefactor_expression() * contracted / denominator
expected = (
    2
    * e**4
    / s**2
    * (
        2 * me**2 * (2 * mm**2 + s - t - u)
        + 2 * me**4
        + 2 * mm**4
        + 2 * mm**2 * (s - t - u)
        + t**2
        + u**2
    )
)
assert (squared - expected).replace(u, 2 * me**2 + 2 * mm**2 - s - t).together() == E(
    "0"
)
massless = squared.replace(me, E("0")).replace(mm, E("0"))
assert (massless - 2 * e**4 * (t**2 + u**2) / s**2).replace(u, -s - t).together() == E(
    "0"
)
print("Generated FeynCalc QED massive and massless squared matrix elements passed")

# Carry the generated matrix element through physical two-body phase space.
# Positive s fixes the square-root branch; Symbolica performs the integral.
s_physical = S("s_physical", is_positive=True)
costheta, alpha = S("cos_theta", "alpha")
pi = Expression.PI
massless_kin = fk.Kinematics.mandelstam(
    [p1, p2, k1, k2], [E("0")] * 4, [s_physical, t, u]
)
flux = massless_kin.flux(p1, p2)
measure = massless_kin.two_body_phase_space(k1, k2)
assert flux == 2 * s_physical
assert (measure - 1 / (32 * pi**2)).together() == E("0")
angular = (
    massless.replace(s, s_physical)
    .replace(t, -s_physical * (1 - costheta) / 2)
    .replace(u, -s_physical * (1 + costheta) / 2)
    .replace(e**4, (4 * pi * alpha) ** 2)
)
differential = (angular * measure / flux).expand()
assert (differential - alpha**2 * (1 + costheta**2) / (4 * s_physical)).together() == E(
    "0"
)
primitive = differential.to_polynomial().integrate(costheta).to_expression()
total = (
    2
    * pi
    * (primitive.replace(costheta, E("1")) - primitive.replace(costheta, E("-1")))
)
assert (total - 4 * pi * alpha**2 / (3 * s_physical)).together() == E("0")

# Independently use t as the polar coordinate: d(cos(theta))/dt = 2/s.
differential_t = (
    2 * differential.replace(costheta, 1 + 2 * t / s_physical) / s_physical
).expand()
expected_t = alpha**2 * (s_physical**2 + 2 * s_physical * t + 2 * t**2) / s_physical**4
assert (differential_t - expected_t).together() == 0
primitive_t = differential_t.to_polynomial().integrate(t).to_expression()
total_t = 2 * pi * (primitive_t.replace(t, 0) - primitive_t.replace(t, -s_physical))
assert (total_t - total).together() == 0

# The same normalization gives the isotropic scalar two-body rest-frame width.
parent, daughter1, daughter2 = S("parent", "daughter1", "daughter2")
mass = S("parent_mass", is_positive=True)
decay_kin = (
    fk.Kinematics()
    .with_scalar_product(parent, parent, mass**2)
    .with_scalar_product(daughter1, daughter1, E("0"))
    .with_scalar_product(daughter2, daughter2, E("0"))
    .with_scalar_product(daughter1, daughter2, mass**2 / 2)
)
width_per_matrix_element = (
    4
    * pi
    * decay_kin.two_body_phase_space(daughter1, daughter2)
    / decay_kin.flux(parent)
)
assert (width_per_matrix_element - 1 / (16 * pi * mass)).together() == E("0")
assert (
    fk.FourMomentum(5.0, 0.0, 0.0, 4.0).flux(fk.FourMomentum(13.0, 0.0, 0.0, -12.0))
    == 448.0
)
assert fk.FourMomentum(13.0, 0.0, 0.0, 12.0).flux() == 26.0
print(
    "Generated QED angular distribution and total cross section, scalar decay normalization passed"
)

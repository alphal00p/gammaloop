"""Generated different-flavor QCD annihilation, including spin and color averages.

Reference: https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QiQibar-QjQjbar
Bottom and top flavors keep both quark masses symbolic in the SM fixture.
"""

from pathlib import Path

from symbolica import E, S, Symbol
from symbolica.community import hep as fk
from symbolica.community.spenso import ColorSimplifySettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
result = (
    fk.Process(model, [5, -5], [6, -6])
    .with_loop_count(1, 1)
    .generate_cross_section(
        max_vertices=4,
        maximum_bridges=None,
        vertex_allow=["V_76", "V_137"],
        numerator_grouping=None,
        progress=None,
    )
)
assert len(result.diagrams) == 1
diagram = result.diagrams[0]
assert len(diagram.cuts) == 1
cut = diagram.cuts[0]
assert sorted(p.pdg_code for p in cut.particles) == [-6, 6]
assert all(side.loop_count == 0 for side in (cut.left, cut.right))
Q, P, K = S("gammalooprs::Q", "gammalooprs::P", "gammalooprs::K")
s, t, u, mb, mt, gs, k1, k2 = S(
    "s", "t", "u", "UFO::MB", "UFO::MT", "UFO::G", "k1", "k2"
)
kinematics = fk.Kinematics.mandelstam(
    [P(0), P(1), k1, k2], [mb**2, mb**2, mt**2, mt**2], [s, t, u]
)
projector = diagram.projector_expression()
for edge in diagram.external_edges:
    projector = model.particle_by_pdg(edge.particle_pdg).sum_spins(
        projector, Q(edge.id), edge=edge.id, average=True
    )
numerator = model.expand_couplings(
    diagram.numerator_expression().to_expression() * projector
)

# Associate open color slots with the generated external wavefunction indices.
# Slot duals determine the closure orientation, including the antiquark; no
# fixed half-edge numbering or hand-written SU(3) averaging factor is needed.
bis, wave, index, metric, left, right = S(
    "spenso::bis", "wave_", "index_", "spenso::g", "left_", "right_"
)
slots = TensorExpression(numerator).structure.slots
color_projector = E("1")
initial_color_states = 1
for edge in diagram.external_edges:
    particle = model.particle_by_pdg(edge.particle_pdg)
    initial_color_states *= abs(particle.color)
    indices = [
        dict(match)[index]
        for match in diagram.projector_expression().match(
            wave(edge.id, bis(4, index)), max_level=0
        )
    ]
    closure_slots = [
        slot.dual().to_expression()
        for slot in slots
        if any(
            slot.to_expression().replace(i, E("0")) != slot.to_expression()
            for i in indices
        )
    ]
    assert len(closure_slots) == 2
    closure = metric(*closure_slots)
    color_indices = dict(
        next(closure.match(particle.color_sum(left, right), max_level=0))
    )
    color_projector *= particle.color_sum(
        color_indices[left], color_indices[right], average=True
    )

# Check singlet and adjoint sums through the same installed public interface.
i, j = S("color_i", "color_j")
assert model.particle_by_pdg(11).color_sum(i, j, average=True) == E("1")
for pdg in (5, -5, 21):
    identity = model.particle_by_pdg(pdg).color_sum(i, i, average=True)
    assert TensorExpression(identity).simplify_metrics().to_expression() == E("1")

# A named adjoint dimension is required by Spenso. Impose dA=Nc²-1 after
# contraction, when the dimensions are ordinary scalar expressions.
Nc, dA, cof, coad = S("Nc", "dA", "spenso::cof", "spenso::coad")
generic = (
    (numerator * color_projector * initial_color_states / Nc**2)
    .replace(cof(3, index), cof(Nc, index))
    .replace(coad(8, index), coad(dA, index))
)
contracted = (
    TensorExpression(generic.expand())
    .simplify_color(ColorSimplifySettings(substitute_cof_dimension_invariants=True))
    .simplify_gamma()
    .expand()
    .to_dots()
)
assert contracted.is_scalar
loop_edge = diagram.loop_momentum_basis.loop_edges[0]
loop_particle = next(
    p.pdg_code for edge, p in zip(cut.edges, cut.particles) if edge.id == loop_edge
)
physical_momentum = k1 if loop_particle == 6 else k2
routed = diagram.loop_momentum_basis.route_expression(
    contracted.to_expression()
).replace(K(0, index), cut.orientations[loop_edge] * physical_momentum(index))
a, b, c, inverse, denom = S("a_", "b_", "c_", "inverse_", "gammalooprs::denom")
denominator = kinematics.apply(
    diagram.denominator_expression(
        edge_powers={edge.id: 0 for edge in cut.edges}, dimension=4, in_lmb=True
    )
    .to_expression()
    .replace(denom(a, b, c, inverse), inverse)
)
assert (denominator - s**2).expand() == E("0")
squared = (
    diagram.overall_factor_expression(evaluate=True)
    * diagram.numerator_prefactor_expression()
    * kinematics.apply(routed)
    / denominator
).replace(dA, Nc**2 - 1)
expected = (
    (Nc**2 - 1)
    * gs**4
    / (2 * Nc**2 * s**2)
    * (
        2 * mb**2 * (2 * mt**2 + s - t - u)
        + 2 * mb**4
        + 2 * mt**4
        + 2 * mt**2 * (s - t - u)
        + t**2
        + u**2
    )
)
assert (squared - expected).replace(u, 2 * mb**2 + 2 * mt**2 - s - t).together() == E(
    "0"
)
massless = squared.replace(mb, E("0")).replace(mt, E("0"))
assert (massless - (Nc**2 - 1) * gs**4 * (t**2 + u**2) / (2 * Nc**2 * s**2)).replace(
    u, -s - t
).together() == E("0")
su3 = massless.replace(Nc, E("3"))
assert (su3 - 4 * gs**4 * (t**2 + u**2) / (9 * s**2)).replace(
    u, -s - t
).together() == E("0")

# Carry the generated massless SU(3) result through physical two-body phase space.
s_physical = S("s_physical", is_positive=True)
costheta, alpha_s = S("costheta", "alpha_s")
pi = Symbol.PI
massless_kinematics = fk.Kinematics.mandelstam(
    [P(0), P(1), k1, k2], [E("0")] * 4, [s_physical, t, u]
)
angular = (
    su3.replace(s, s_physical)
    .replace(t, -s_physical * (1 - costheta) / 2)
    .replace(u, -s_physical * (1 + costheta) / 2)
    .replace(gs**4, (4 * pi * alpha_s) ** 2)
)
differential = (
    angular
    * massless_kinematics.two_body_phase_space(k1, k2)
    / massless_kinematics.flux(P(0), P(1))
).expand()
assert (
    differential - alpha_s**2 * (1 + costheta**2) / (18 * s_physical)
).together() == E("0")
primitive = differential.to_polynomial().integrate(costheta).to_expression()
total = (
    2
    * pi
    * (primitive.replace(costheta, E("1")) - primitive.replace(costheta, E("-1")))
)
assert (total - 8 * pi * alpha_s**2 / (27 * s_physical)).together() == E("0")
print(
    "Generated QCD annihilation: massive SU(N), massless SU(3), angular distribution and total cross section passed"
)

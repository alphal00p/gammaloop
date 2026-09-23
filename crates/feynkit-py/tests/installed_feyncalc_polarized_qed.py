"""Generated chiral-projected QED production, independent of cut momentum choice.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-MuAmu
Chirality coincides with the specified helicities only in the massless limit.
"""

from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import Representation, TensorExpression, chain

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
process = fk.Process.cross_section([11, -11], [13, -13]).with_loop_count(1, 1)
result = fk.Generator(model).generate(
    process,
    max_vertices=4,
    maximum_bridges=None,
    vertex_allow=["V_98", "V_99"],
    numerator_grouping=None,
    progress=None,
)
assert len(result.diagrams) == 1
generated = result.diagrams[0]
Q, P, K = S("gammalooprs::Q", "gammalooprs::P", "gammalooprs::K")
s, t, u, me, mm, e = S("s", "t", "u", "UFO::Me", "UFO::MM", "UFO::ee")
a, b, ell = S("chiral::a_", "chiral::b_", "chiral::ell_")
gamma = TensorExpression.gamma(4)(a, b, ell).to_expression()
right_current = chain(
    Representation.bis(4)(a),
    Representation.bis(4)(b),
    TensorExpression.gamma(4)(a, "middle", ell),
    TensorExpression.projp(4)("middle", b),
).to_expression()

# Both physical final-state particles can serve as the loop coordinate.
# Its charge comes from the cut, not from the stored edge's particle record.
for coordinate in generated.cuts[0].edges:
    diagram = generated.with_loop_momentum_edges([coordinate.id])
    cut = diagram.cuts[0]
    particles = {
        edge.id: particle.pdg_code for edge, particle in zip(cut.edges, cut.particles)
    }
    assert sorted(particles.values()) == [-13, 13]
    orientation = cut.orientations[coordinate.id]
    physical, other, index = S(
        "selected_final_momentum", "other_final_momentum", "index_"
    )
    k1, k2 = (physical, other) if particles[coordinate.id] == 13 else (other, physical)
    kin = fk.Kinematics.mandelstam(
        [P(0), P(1), k1, k2], [me**2, me**2, mm**2, mm**2], [s, t, u]
    )

    numerator = E("1")
    for vertex in diagram.vertices:
        numerator *= (
            vertex.numerator_expression().to_expression().replace(gamma, right_current)
        )
    for edge in diagram.internal_edges:
        numerator *= edge.numerator_expression().to_expression()
    projector = diagram.projector_expression()
    for edge in diagram.external_edges:
        # Specified initial states carry no spin average.
        projector = model.particle_by_pdg(edge.particle_pdg).sum_spins(
            projector, Q(edge.id), edge=edge.id
        )
    numerator = model.expand_couplings(numerator * projector)
    contracted = (
        TensorExpression(numerator.expand())
        .simplify_gamma()
        .expand()
        .simplify_epsilon()
        .expand()
        .to_dots()
        .to_expression()
    )
    contracted = kin.apply(
        diagram.loop_momentum_basis.route_expression(contracted).replace(
            K(0, index), orientation * physical(index)
        )
    )
    denominator = diagram.denominator_expression(
        edge_powers={edge.id: 0 for edge in cut.edges}, dimension=4, in_lmb=True
    ).to_expression()
    q, m, power, inverse = S("den::q_", "den::m_", "den::power_", "den::inverse_")
    denominator = kin.apply(
        denominator.replace(
            S("gammalooprs::denom")(q, m, power, inverse), inverse
        ).replace(K(0, index), orientation * physical(index))
    )
    assert (denominator - s**2).together() == E("0")
    squared = (
        diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
        * contracted
        / denominator
    ).together()
    expected = 4 * e**4 * (me**2 + mm**2 - u) ** 2 / s**2
    assert (squared - expected).replace(
        u, 2 * me**2 + 2 * mm**2 - s - t
    ).together() == E("0")
    massless = squared.replace(me, E("0")).replace(mm, E("0"))
    assert (massless - 4 * e**4 * u**2 / s**2).replace(u, -s - t).together() == E("0")
    print(
        f"Generated polarized QED passed with PDG {particles[coordinate.id]} as loop coordinate"
    )

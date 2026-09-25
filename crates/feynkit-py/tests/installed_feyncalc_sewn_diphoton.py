"""Massive diphoton annihilation from physical cuts of sewn forward graphs.

Compare both cut-momentum coordinates, three photon projectors and Ward identities
with the ordinary-amplitude reference in installed_feyncalc_diphoton.py.
Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-GaGa
"""

from symbolica import E, Replacement, S
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
electron, photon = (model.particle(name) for name in ("e-", "a"))
vertices = [
    vertex
    for vertex in model.vertex_rules
    if sorted(vertex.particles)
    == sorted([electron.name, electron.antiname, photon.name])
]
assert len(vertices) == 1
result = model.process(
    ["e-", "e+"], ["a", "a"], vertex_allow=vertices
).generate_cross_section(
    loops=1,
    max_vertices=4,
    maximum_bridges=None,
    numerator_grouping=None,
    progress=None,
)
assert len(result.diagrams) == 2
Q, P, K = S("gammalooprs::Q", "gammalooprs::P", "gammalooprs::K")
s, t, u, mass, charge = S("s", "t", "u", "UFO::Me", "UFO::ee")
k1, k2, index = S("k1", "k2", "index_")
a, b, c, inverse = S("a_", "b_", "c_", "inverse_")
metric, mink = S("spenso::g", "spenso::mink")
kin = hep.Kinematics.mandelstam(
    [P(0), P(1), k1, k2], [mass**2, mass**2, E("0"), E("0")], [s, t, u]
)
x, y = t - mass**2, u - mass**2
expected = (
    2
    * charge**4
    * (x / y + y / x + 4 * mass**2 * s / (x * y) - 4 * mass**4 * s**2 / (x**2 * y**2))
)
expected = expected.replace(s, 2 * mass**2 - t - u).together()
results = {}
for coordinate_index in (0, 1):
    for mode in ("covariant", "null", "timelike", "Ward 0", "Ward 1"):
        total = E("0")
        for original in result.diagrams:
            coordinate = original.cuts[0].edges[coordinate_index].id
            diagram = original.with_loop_momentum_edges([coordinate])
            cut = diagram.cuts[0]
            assert [particle.name for particle in cut.particles] == ["a", "a"]
            assert all(side.loop_count == 0 for side in (cut.left, cut.right))
            orientation = cut.orientations[coordinate]
            projector = diagram.projector_expression()
            for edge in diagram.external_edges:
                projector = model.particle(edge.particle_name).sum_spins(
                    projector, Q(edge.id), edge=edge.id, average=True
                )
            numerator = model.expand_couplings(
                diagram.numerator_expression().to_expression() * projector
            )
            for position, edge in enumerate(cut.edges):
                if mode == "covariant" or (
                    mode.startswith("Ward") and mode != f"Ward {position}"
                ):
                    continue
                # The photon propagator numerator is -i*g. Replacing g by the
                # negative physical density retains its original propagator phase.
                edge_metric = E("1𝑖") * edge.numerator_expression().to_expression()
                match = dict(
                    next(edge_metric.match(metric(mink(4, a), mink(4, b)), max_level=0))
                )
                if mode.startswith("Ward"):
                    physical = Q(edge.id, mink(4, match[a])) * Q(
                        edge.id, mink(4, match[b])
                    )
                else:
                    reference = (
                        P(0) if mode == "timelike" else Q(cut.edges[1 - position].id)
                    )
                    physical = photon.spin_sum(
                        Q(edge.id), match[a], match[b], reference=reference
                    )
                numerator = numerator.replace(edge_metric, -physical)
            scalar = (
                TensorExpression(
                    numerator.expand(), cook_indices=CookSettings.indices()
                )
                .simplify_gamma()
                .expand()
                .simplify_gamma()
                .expand()
                .simplify_metrics()
                .to_dots()
            )
            assert scalar.is_scalar
            # A cut coordinate is a stored directed edge momentum. Convert it
            # to the positive-energy outgoing photon before imposing kinematics.
            contracted = kin.apply(
                diagram.loop_momentum_basis.route_expression(
                    scalar.to_expression()
                ).replace(K(0, index), orientation * k1(index))
            )
            denominator = kin.apply(
                diagram.denominator_expression(
                    edge_powers={edge.id: 0 for edge in cut.edges},
                    dimension=4,
                    in_lmb=True,
                )
                .to_expression()
                .replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
                .replace(K(0, index), orientation * k1(index))
            )
            total += (
                diagram.overall_factor_expression(evaluate=True)
                * diagram.numerator_prefactor_expression()
                * contracted
                / denominator
            ).replace(s, 2 * mass**2 - t - u)
        # The forward graphs identify the two identical cut photons. Restore
        # both labeled final-state assignments for the labeled squared amplitude.
        squared = (
            total + total.replace_multiple([Replacement(t, u), Replacement(u, t)])
        ).together()
        target = E("0") if mode.startswith("Ward") else expected
        assert (squared - target).together() == 0, (coordinate_index, mode)
        if not mode.startswith("Ward"):
            assert (
                squared.replace(mass, 0) - 2 * charge**4 * (t / u + u / t)
            ).together() == 0
            assert squared.replace(mass, 1).replace(charge, 1).replace(t, -1).replace(
                u, -7
            ).together() == E("83/8")
        results[coordinate_index, mode] = squared
        print("PASS", coordinate_index, mode, flush=True)
assert results[0, "covariant"] == results[1, "covariant"]
print(
    "Both cut coordinates, physical photon sums, Ward identities and massive diphoton reference passed",
    flush=True,
)

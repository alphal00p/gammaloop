"""Massive sewn QED benchmarks with physical spin sums and cut routing.

References: FeynCalc QED/Tree ElGa-ElGa, ElAel-MuAmu, and ElMu-ElMu.
The Compton cases check both covariant and timelike axial incoming-photon sums.
"""

from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
Q, P, K = S("gammalooprs::Q", "gammalooprs::P", "gammalooprs::K")
s, t, u, me, mm, e = S("s", "t", "u", "UFO::Me", "UFO::MM", "UFO::ee")
k1, k2, index = S("k1", "k2", "index_")
a, b, c, inverse = S("a_", "b_", "c_", "inverse_")
electron, muon, photon = (model.particle_by_pdg(pdg) for pdg in (11, 13, 22))
vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles)
    in (
        sorted([electron.antiname, electron.name, photon.name]),
        sorted([muon.antiname, muon.name, photon.name]),
    )
]
assert len(vertices) == 2

for name, incoming, outgoing in [
    ("electron Compton", [11, 22], [11, 22]),
    ("positron Compton", [-11, 22], [-11, 22]),
    ("annihilation", [11, -11], [13, -13]),
    ("electron-muon scattering", [11, 13], [11, 13]),
]:
    compton = name.endswith("Compton")
    masses = (
        [me**2, E("0"), me**2, E("0")]
        if compton
        else [me**2, me**2, mm**2, mm**2]
        if name == "annihilation"
        else [me**2, mm**2, me**2, mm**2]
    )
    kin = fk.Kinematics.mandelstam([P(0), P(1), k1, k2], masses, [s, t, u])
    result = fk.Process(model, incoming, outgoing).generate_cross_section(
        loops=1,
        max_vertices=4,
        maximum_bridges=None,
        vertex_allow=vertices,
        numerator_grouping=None,
        progress=None,
    )
    assert len(result.diagrams) == (4 if compton else 1)
    if compton:
        expected = (
            2
            * e**4
            * (
                -(me**4) * (3 * s**2 + 14 * s * u + 3 * u**2)
                + me**2 * (s**3 + 7 * s**2 * u + 7 * s * u**2 + u**3)
                + 6 * me**8
                - s * u * (s**2 + u**2)
            )
            / ((s - me**2) ** 2 * (u - me**2) ** 2)
        )
    elif name == "annihilation":
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
    else:
        expected = (
            2
            * e**4
            / t**2
            * (
                -2 * me**2 * (-2 * mm**2 + s - t + u)
                + 2 * me**4
                + 2 * mm**4
                - 2 * mm**2 * (s - t + u)
                + s**2
                + u**2
            )
        )
    for reference in [None, P(0)] if compton else [None]:
        squared = E("0")
        for diagram in result.diagrams:
            assert len(diagram.cuts) == 1
            cut = diagram.cuts[0]
            assert sorted(p.pdg_code for p in cut.particles) == sorted(outgoing)
            coordinate = next(
                edge.id
                for edge, particle in zip(cut.edges, cut.particles)
                if particle.pdg_code == outgoing[0]
            )
            diagram = diagram.with_loop_momentum_edges([coordinate])
            cut = diagram.cuts[0]
            assert all(side.loop_count == 0 for side in (cut.left, cut.right))
            orientation = cut.orientations[coordinate]
            projector = diagram.projector_expression()
            for edge in diagram.external_edges:
                projector = model.particle_by_pdg(edge.particle_pdg).sum_spins(
                    projector,
                    Q(edge.id),
                    edge=edge.id,
                    average=True,
                    reference=reference,
                    covariant=reference is None,
                )
            numerator = model.expand_couplings(
                diagram.numerator_expression().to_expression() * projector
            )
            numerator = (
                TensorExpression(numerator.expand())
                .simplify_gamma()
                .expand()
                .to_dots()
                .to_expression()
            )
            numerator = kin.apply(
                diagram.loop_momentum_basis.route_expression(numerator).replace(
                    K(0, index), orientation * k1(index)
                )
            )
            denominator = diagram.denominator_expression(
                edge_powers={edge.id: 0 for edge in cut.edges},
                dimension=4,
                in_lmb=True,
            ).to_expression()
            denominator = kin.apply(
                denominator.replace(
                    S("gammalooprs::denom")(a, b, c, inverse), inverse
                ).replace(K(0, index), orientation * k1(index))
            )
            squared += (
                diagram.overall_factor_expression(evaluate=True)
                * diagram.numerator_prefactor_expression()
                * numerator
                / denominator
            )
        difference = (squared - expected).replace(t, sum(masses) - s - u).together()
        assert difference == 0, (name, reference, difference)
        if compton:
            massless = squared.replace(t, 2 * me**2 - s - u).replace(me, E("0"))
            assert (massless + 2 * e**4 * (s / u + u / s)).together() == 0
        print(
            f"Sewn {name}, reference={reference}: exact massive result passed",
            flush=True,
        )

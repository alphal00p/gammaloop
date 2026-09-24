"""Generated Compton amplitudes, their adjoints, and physical spin sums.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElGa-ElGa
The sewn-forward-graph conversion is a separate regression requirement.
"""

from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
s, t, u, mass, charge = S("s", "t", "u", "UFO::Me", "UFO::ee")
kinematics = fk.Kinematics.mandelstam(
    [P(0), P(1), P(2), P(3)], [mass**2, E("0"), mass**2, E("0")], [s, t, u]
)
a, b, c, inverse = S("a_", "b_", "c_", "inverse_")
wave, representation, index = S("wave_", "representation_", "index_")
ports = S("external_0", "external_1", "external_2", "external_3")
conjugate, adjoint_index = S("spenso::conj", "adjoint_index")
expected = (
    2
    * charge**4
    * (
        -(mass**4) * (3 * s**2 + 14 * s * u + 3 * u**2)
        + mass**2 * (s**3 + 7 * s**2 * u + 7 * s * u**2 + u**3)
        + 6 * mass**8
        - s * u * (s**2 + u**2)
    )
    / ((s - mass**2) ** 2 * (u - mass**2) ** 2)
)

for pdg in (11, -11):
    generated = fk.Generator(model).generate(
        fk.Process.amplitude([pdg, 22], [pdg, 22]).with_loop_count(0, 0),
        max_vertices=2,
        maximum_bridges=None,
        vertex_allow=["V_98"],
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 2
    amplitude = E("0")
    for diagram in generated.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
        # Align external ports using the generated wavefunctions, independently
        # of the internal half-edge numbering in the two diagrams.
        for edge in diagram.external_edges:
            matches = list(
                diagram.projector_expression().match(
                    wave(edge.id, representation(4, index)), max_level=0
                )
            )
            assert len(matches) == 1
            matched = dict(matches[0])
            numerator = numerator.replace(
                matched[representation](4, matched[index]),
                matched[representation](4, ports[edge.external_index]),
            )
        denominator = diagram.denominator_expression(
            dimension=4, in_lmb=True
        ).to_expression()
        denominator = kinematics.apply(
            denominator.replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
        )
        amplitude += (
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
            / denominator
        )

    operator = TensorExpression(amplitude.expand())
    assert len(operator.structure.slots) == 4
    adjoint = operator.dirac_adjoint().expand().simplify_gamma0().to_expression()
    # Physical external momenta, the tree-level charge, masses and invariants
    # are real. Idenso retains these assumptions as explicit conjugations.
    adjoint = adjoint.replace(conjugate(P(a, b)), P(a, b))
    for real in (mass, charge, s, t, u):
        adjoint = adjoint.replace(conjugate(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(adjoint_index)

    for axial in (False, True):
        density = E("1")
        # Dirac adjunction exchanges matrix input/output roles. Completeness
        # tensors join the original ket to the corresponding adjoint bra.
        for particle, momentum, left, right, average in [
            (pdg, P(0), ports[0], adjoint_index(ports[2]), True),
            (pdg, P(2), adjoint_index(ports[0]), ports[2], False),
            (22, P(1), ports[1], adjoint_index(ports[1]), True),
            (22, P(3), ports[3], adjoint_index(ports[3]), False),
        ]:
            if particle < 0:
                left, right = right, left
            density *= model.particle_by_pdg(particle).spin_sum(
                momentum,
                left,
                right,
                average=average,
                reference=P(0) if axial and particle == 22 else None,
            )
        squared = (
            TensorExpression(
                operator.to_expression() * adjoint * density,
                cook_indices=CookSettings.indices(),
            )
            .expand()
            .simplify_gamma()
            .expand()
            .to_dots()
            .to_expression()
        )
        squared = kinematics.apply(squared).replace(t, 2 * mass**2 - s - u).together()
        assert (squared - expected).together() == E("0")
        assert (
            squared.replace(mass, E("0")) + 2 * charge**4 * (s / u + u / s)
        ).together() == E("0")
        point = squared
        for parameter, value in ((charge, 1), (mass, 1), (s, 3), (u, 0)):
            point = point.replace(parameter, E(str(value)))
        assert point.together() == E("3")
        print(
            f"Generated Compton: PDG={pdg}, axial={axial}, massive and massless references passed"
        )

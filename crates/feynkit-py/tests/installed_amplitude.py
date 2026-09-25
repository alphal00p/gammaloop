"""Coherent amplitudes close physical ports and reproduce a massive QED result."""

from symbolica import E, S
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model.standard_model()
generated = model.process(
    ["e-", "e+"], ["a", "a"], vertex_allow=["V_98"]
).generate_diagrams(
    max_vertices=2, maximum_bridges=None, numerator_grouping=None, progress=None
)
amplitude = hep.Amplitude(generated.diagrams)
assert len(amplitude.diagrams) == len(amplitude.terms) == 2
assert [leg.index for leg in amplitude.legs] == [0, 1, 2, 3]
assert [leg.particle.name for leg in amplitude.legs] == ["e-", "e+", "a", "a"]
assert [leg.state for leg in amplitude.legs] == ["incoming"] * 2 + ["outgoing"] * 2
assert isinstance(amplitude.expression(), TensorExpression)
assert amplitude.expression().rank == 4
assert amplitude.conjugate().is_conjugated
assert not amplitude.conjugate().conjugate().is_conjugated
assert "conj" not in str(amplitude.conjugate().expression().to_expression())
assert amplitude._repr_html_()

squared = amplitude.squared()
assert squared.expression().rank == 8
assert squared.spin_summed == []
fermions = squared.sum_spins([0, 1], average_initial=True)
assert fermions.expression().rank == 4
assert fermions.spin_summed == [0, 1]

P = S("gammalooprs::P")
mass, charge = S("UFO::Me", "UFO::ee")
s, t, u = S("amplitude_test::s", "amplitude_test::t", "amplitude_test::u")
kin = hep.Kinematics.mandelstam(
    [P(i) for i in range(4)], [mass**2, mass**2, E("0"), E("0")], [s, t, u]
)
x, y = t - mass**2, u - mass**2
expected = (
    2
    * charge**4
    * (x / y + y / x + 4 * mass**2 * s / (x * y) - 4 * mass**4 * s**2 / (x**2 * y**2))
)
for references in (None, {2: P(3), 3: P(2)}, {2: P(0), 3: P(0)}):
    closed = fermions.sum_spins(references=references).sum_colors()
    expression = closed.expression()
    assert expression.is_scalar
    scalar = (
        expression.expand()
        .simplify_gamma()
        .expand()
        .simplify_gamma()
        .simplify_metrics()
        .to_dots()
    )
    result = kin.apply(scalar.to_expression())
    assert (result - expected).replace(s, 2 * mass**2 - t - u).together() == E("0")
    assert closed.sum_spins().expression().to_expression() == expression.to_expression()

for operation in (
    lambda: hep.Amplitude([]),
    lambda: fermions.sum_spins([0], average_initial=True),
    lambda: squared.sum_colors([99]),
    lambda: squared.sum_spins([0], references={3: P(1)}),
):
    try:
        operation()
    except hep.AmplitudeError:
        pass
    else:
        raise AssertionError("invalid amplitude operation was accepted")

assert hep.SquaredAmplitude.from_diagram(generated.diagrams[0]).expression().rank == 8
D = S("amplitude_test::D")
dimensional = hep.Amplitude(generated.diagrams, dimension=D)
assert [leg.slots[0].representation.dimension for leg in dimensional.legs] == [
    E("4"),
    E("4"),
    D,
    D,
]
assert dimensional.squared().sum_spins().expression().is_scalar

z = S("amplitude_test::z")
scalar_diagram = hep.FeynmanDiagram.from_dot(
    model,
    'digraph { a [num="amplitude_test::z"]; ext [style=invis]; '
    'ext -> a [particle="H"]; a -> ext [particle="H"]; a -> ext [particle="H"]; }',
)
real_amplitude = hep.Amplitude.from_diagram(scalar_diagram, real=[z])
assert real_amplitude.conjugate().expression().to_expression() == z
assert real_amplitude.squared().expression().to_expression() == z**2

print(
    "Amplitude: physical ports, conjugation, interference, state sums and massive diphoton oracle passed"
)

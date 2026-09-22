"""Generated QED Pauli form factor, shared IBP and the Schwinger limit.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/El-GaEl
The on-shell electron vertex is projected at generic spacelike photon momentum
before taking t->0, where the two form-factor projectors become degenerate.
"""

from pathlib import Path

import numpy as np
from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
dimension, transfer, mass, charge, mu = S(
    "gminus2::D", "gminus2::t", "UFO::Me", "UFO::ee", "gminus2::mu"
)
K, P, mink, bis, metric, gamma = S(
    "gammalooprs::K",
    "gammalooprs::P",
    "spenso::mink",
    "spenso::bis",
    "spenso::g",
    "spenso::gamma",
)
a, b, c, d, index, wave, slot = S(
    "gminus2::a", "gminus2::b", "gminus2::c", "gminus2::d", "index_", "wave_", "slot_"
)
kinematics = (
    hep.Kinematics(dimension, momenta=[K(0), P(0), P(1)])
    .with_scalar_product(P(0), P(0), mass**2)
    .with_scalar_product(P(1), P(1), transfer)
    .with_scalar_product(P(0), P(1), transfer / 2)
)
operators, diagrams = [], []
for loops in (0, 1):
    generated = hep.Generator(model).generate(
        hep.Process.amplitude([11], [22, 11]).with_loop_count(loops, loops),
        max_vertices=1 + 2 * loops,
        maximum_bridges=0,
        vertex_allow=["V_98"],
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    operator = (
        model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    for edge in diagram.external_edges:
        representation = mink if edge.external_index == 1 else bis
        matches = list(
            diagram.projector_expression().match(
                wave(edge.id, representation(4, index)), max_level=0
            )
        )
        assert len(matches) == 1
        port = dict(matches[0])[index]
        operator = operator.replace(
            representation(4, port),
            representation(4, {0: b, 1: mu, 2: a}[edge.external_index]),
        )
    operators.append(operator.replace(mink(4, index), mink(dimension, index)))
    diagrams.append(diagram)

# Normalize to the generated tree vertex, including its external fermion sign.
# A scalar loop contributes i/(16*pi^2) in the OneLOop convention below.
tree_coupling = (
    operators[0] / gamma(bis(4, a), bis(4, b), mink(dimension, mu))
).together()
assert (tree_coupling - Symbol.I * charge).expand() == E("0")
vertex = (Symbol.I * operators[1] / (tree_coupling * charge**2)).expand()
particle = model.particle_by_pdg(11)
incoming = particle.spin_sum(P(0), b, c)
outgoing = particle.spin_sum(P(2), d, a).replace(P(2, slot), P(0, slot) - P(1, slot))
incoming = incoming.replace(mink(4, index), mink(dimension, index)).replace(
    mink(4), mink(dimension)
)
outgoing = outgoing.replace(mink(4, index), mink(dimension, index)).replace(
    mink(4), mink(dimension)
)
vector_sum = 2 * P(0, mink(dimension, mu)) - P(1, mink(dimension, mu))
basis = [
    gamma(bis(4, a), bis(4, b), mink(dimension, mu)),
    vector_sum * metric(bis(4, a), bis(4, b)) / (2 * mass),
]
probes = [
    gamma(bis(4, c), bis(4, d), mink(dimension, mu)),
    vector_sum * metric(bis(4, c), bis(4, d)) / (2 * mass),
]
traces = []
for operator in basis + [vertex]:
    for probe in probes:
        trace = (
            TensorExpression((operator * incoming * probe * outgoing).expand())
            .simplify_gamma()
            .expand()
            .to_dots()
            .to_expression()
        )
        traces.append(kinematics.apply(trace).expand())
# These independently specified traces also fix the normalization of the
# projector basis gamma^mu and (p_in+p_out)^mu/(2m).
expected_gram = [
    2 * (dimension - 2) * transfer + 8 * mass**2,
    8 * mass**2 - 2 * transfer,
    8 * mass**2 - 2 * transfer,
    (4 * mass**2 - transfer) ** 2 / (2 * mass**2),
]
assert all(
    (actual - expected).expand() == E("0")
    for actual, expected in zip(traces[:4], expected_gram, strict=True)
)
projectors = Matrix.from_linear(2, 2, traces[:4])
pauli_numerator = projectors.solve(Matrix.vec(traces[4:]))[1, 0].to_expression()
family = diagrams[1].integral_family(kinematics=kinematics)
denominators = S("gminus2::d0", "gminus2::d1", "gminus2::d2")
polynomial = family.rewrite_numerator(pauli_numerator, denominators).expand()
terms = []
for monomial, coefficient in polynomial.coefficient_list(*denominators):
    powers = [
        1 - monomial.to_polynomial(vars=denominators).degree(label)
        for label in denominators
    ]
    terms.append((powers, coefficient))
solution = hep.IBPFamily(family, name="electron_vertex").reduce_laporta(
    [powers for powers, coefficient in terms], max_depth=2
)
assert {tuple(powers) for powers in solution.residuals} == {
    (0, 1, 0),
    (1, 0, 0),
    (1, 1, 0),
}
integral, tadpole, bubble = S("gminus2::I", "gminus2::A0", "gminus2::B0")
reduction = sum(
    (
        coefficient * solution.reduce(powers, integral=integral)
        for powers, coefficient in terms
    ),
    E("0"),
)
# Both one-line pinches are massive tadpoles; their loop shifts have unit Jacobian.
# The two-line pinch carries q^2=t and equal masses m^2.
reduction = reduction.replace_multiple(
    [
        Replacement(integral(0, 1, 0), tadpole),
        Replacement(integral(1, 0, 0), tadpole),
        Replacement(integral(1, 1, 0), bubble),
    ]
).together()
expected = (
    2
    * (dimension - 5)
    * (-(dimension - 2) * tadpole + 2 * mass**2 * (dimension - 3) * bubble)
    / ((dimension - 3) * (transfer - 4 * mass**2))
)
assert (reduction - expected).together() == E("0")
epsilon, finite_a, finite_b, logarithm = S(
    "gminus2::eps", "gminus2::Af", "gminus2::Bf", "gminus2::L"
)
laurent = (
    reduction.replace(dimension, 4 - 2 * epsilon)
    .replace_multiple(
        [
            Replacement(tadpole, finite_a + mass**2 / epsilon),
            Replacement(bubble, finite_b + 1 / epsilon),
        ]
    )
    .series(epsilon, 0, 0)
    .to_expression()
    .expand()
)
assert laurent.coefficient(epsilon**-1).together() == E("0")
finite = dict(laurent.coefficient_list(epsilon))[E("1")].together()
assert (
    finite - 4 * (finite_a + mass**2 - mass**2 * finite_b) / (transfer - 4 * mass**2)
).together() == E("0")
# Gordon decomposition gives F2=-b. Restoring e^2/(16*pi^2), the b=-2
# limit gives F2(0)=e^2/(8*pi^2)=alpha/(2*pi), independently of m and mu.
limit = (
    finite.replace(transfer, E("0"))
    .replace_multiple(
        [
            Replacement(finite_a, mass**2 * (1 - logarithm)),
            Replacement(finite_b, -logarithm),
        ]
    )
    .together()
)
assert limit == E("-2")
print(
    "Generated vertex: exact Pauli projection, IBP reduction and Schwinger limit passed",
    flush=True,
)

nodes, weights = np.polynomial.legendre.leggauss(96)
nodes, weights = (nodes + 1) / 2, weights / 2
for invariant, electron_mass, scale in [
    (-0.1, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
    (-10.0, 1.0, 1.0),
    (-100.0, 1.0, 1.0),
    (-3.0, 2.0, 7.0),
]:
    # F2(t)/(alpha/(2*pi)); the parameter formula has no IBP/OneLOop input.
    normalized = (
        -complex(
            finite.evaluate(
                {
                    transfer: invariant,
                    mass: electron_mass,
                    finite_a: complex(oneloop.A0(electron_mass**2, scale)[0]),
                    finite_b: complex(
                        oneloop.B0(
                            invariant, electron_mass**2, electron_mass**2, scale
                        )[0]
                    ),
                }
            )
        )
        / 2
    )
    parameter_integral = float(
        np.dot(
            weights,
            electron_mass**2 / (electron_mass**2 - invariant * nodes * (1 - nodes)),
        )
    )
    derivative = 1 + invariant * complex(
        oneloop.dB0(invariant, electron_mass**2, electron_mass**2, scale)[0]
    )
    assert abs(normalized - parameter_integral) < 2e-12
    assert abs(normalized - derivative) < 2e-12
print(
    "Five spacelike Pauli form factors agree with quadrature and the OneLOop bubble derivative",
    flush=True,
)

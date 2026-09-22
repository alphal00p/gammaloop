"""Generated one-loop QED MSbar counterterms in a symbolic covariant gauge.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization
The graph UV expansion, traces, tensor projection and native IBP determine the
bare poles. The tadpole pole and local counterterm operators are stated inputs;
this does not generate counterterm diagrams or a subtraction forest.
"""

import json
from pathlib import Path

from symbolica import E, Matrix, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
D, epsilon, M, mUV, coordinate, integral, xi, s, Nf = S(
    "qed_ren::D",
    "qed_ren::eps",
    "qed_ren::M",
    "qed_ren::mUV",
    "qed_ren::d0",
    "qed_ren::I",
    "qed_ren::xi",
    "qed_ren::s",
    "qed_ren::Nf",
)
specification = json.loads(model.to_json())
for propagator in specification["propagators"]:
    if propagator["particle"] == "a":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))-(1-qed_ren::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(specification))
K, P, mass, charge = S("gammalooprs::K", "gammalooprs::P", "UFO::Me", "UFO::ee")
mink, bis, gamma, metric = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "spenso::g"
)
index, dim, wave, mu, nu = S(
    "qed_ren::index_", "qed_ren::dim_", "qed_ren::wave_", "qed_ren::mu", "qed_ren::nu"
)
den, edge_, mom_, mass_, quad_ = S(
    "gammalooprs::denom",
    "qed_ren::edge_",
    "qed_ren::mom_",
    "qed_ren::mass_",
    "qed_ren::quad_",
)
ordering, value = S(
    "feynkit_generator_factor::ExternalFermionOrderingSign", "qed_ren::value_"
)
kinematics = hep.Kinematics(D, momenta=[K(0), P(0)]).with_scalar_product(P(0), P(0), s)
vacuum = hep.Kinematics(D, momenta=[K(0)])
family = hep.IntegralFamily(
    [K(0)], [], [vacuum.scalar_product(K(0), K(0)) - M], kinematics=vacuum
)
reducer = hep.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
gmunu = metric(mink(D, mu), mink(D, nu))
ppmunu = P(0, mink(D, mu)) * P(0, mink(D, nu))
zero, one = E("0"), E("1")
parts, diagrams = {}, {}
for kind, incoming, outgoing, loops in [
    ("tree", [11], [22, 11], 0),
    ("electron", [11], [11], 1),
    ("photon", [22], [22], 1),
    ("vertex", [11], [22, 11], 1),
]:
    generated = hep.Generator(model).generate(
        hep.Process.amplitude(incoming, outgoing).with_loop_count(loops, loops),
        max_vertices=len(incoming) + len(outgoing) - 2 + 2 * loops,
        maximum_bridges=0,
        vertex_allow=["V_98"],
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    diagrams[kind] = diagram
    ports = {}
    for edge in diagram.external_edges:
        rep = (
            mink
            if kind == "photon"
            or (kind in ("tree", "vertex") and edge.external_index == 1)
            else bis
        )
        ports[edge.external_index] = dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, rep(4, index)), max_level=0
                )
            )
        )[index]
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    if loops:
        numerator = diagram.momentum_basis().route_expression(
            diagram.uv_expansion(mUV, numerator=numerator).to_expression()
        )
    numerator = (
        numerator.replace(mink(dim, index), mink(D, index))
        .replace(mink(dim), mink(D))
        .replace(mUV**2, M)
    )
    if kind == "photon":
        numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
            mink(D, ports[1]), mink(D, nu)
        )
    pattern = den(edge_, mom_, mass_, quad_)
    for match in list(numerator.match(pattern)):
        values = dict(match)
        formal = family.rewrite_numerator(values[quad_], [coordinate])
        assert formal == coordinate
        numerator = numerator.replace(
            den(values[edge_], values[mom_], values[mass_], values[quad_]), formal
        )
    factor = (
        diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    if kind == "electron":
        raw = diagram.overall_factor_expression()
        removed = (raw / raw.replace(ordering(value), one)).replace(
            ordering(value), value
        )
        assert removed == -one
        factor /= removed
        probes = [
            (
                "electron_p",
                gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                * P(0, mink(D, mu))
                / (4 * s),
            ),
            ("electron_m", metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass)),
        ]
    elif kind in ("vertex", "tree"):
        probes = [
            (
                kind,
                gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1])) / (4 * D),
            )
        ]
    else:
        probes = [(kind, one)]
    for label, projector in probes:
        trace = (
            TensorExpression((numerator * projector).expand())
            .simplify_gamma()
            .expand()
            .to_dots()
            .to_expression()
        )
        scalar = family.rewrite_numerator(
            kinematics.apply(reducer.reduce(trace)), [coordinate]
        )
        scalar = (scalar * factor).together().expand()
        terms = []
        for monomial, coefficient in scalar.coefficient_list(coordinate):
            power = -int(
                (monomial.derivative(coordinate) * coordinate / monomial).together()
            )
            assert monomial == coordinate ** (-power)
            assert not coefficient.matches(K(S("qed_ren::args___")))
            terms.append(([power], coefficient))
        parts[label] = terms

tree = parts.pop("tree")[0][1]
assert tree == Symbol.I * charge
solution = hep.IBPFamily(family, name="qed_one_loop").reduce_laporta(
    sorted({tuple(p) for terms in parts.values() for p, c in terms}), max_depth=2
)
assert solution.residuals == [[1]]
assert abs(complex(oneloop.A0(1.0, 1.0)[1]) - 1) < 1e-12
reduced, uv_poles = {}, {}
for label, terms in parts.items():
    expression = sum(
        (c * solution.reduce(p, integral=integral) for p, c in terms), zero
    ).together()
    reduced[label] = expression
    normalized = (
        expression * (Symbol.I / tree if label == "vertex" else one) / charge**2
    )
    if label == "photon":
        normalized *= Nf
    pole = (
        normalized.replace(integral(1), M / epsilon)
        .replace(D, 4 - 2 * epsilon)
        .series(epsilon, 0, -1)
        .to_expression()
        .expand()
    )
    pole = pole.replace(mink(4, index), mink(D, index))
    uv_poles[label] = pole
    assert pole.derivative(M).expand() == zero
    assert pole.derivative(mass).expand() == zero
    assert pole.coefficient(epsilon**-2) == zero
# Compare the four UV structures with the published symbolic-gauge result.
expected = [
    xi / epsilon,
    -(xi + 3) / epsilon,
    -4 * Nf * (s * gmunu - ppmunu) / (3 * epsilon),
    xi / epsilon,
]
for label, reference in zip(
    ("electron_p", "electron_m", "photon", "vertex"), expected, strict=True
):
    assert (uv_poles[label] - reference).together() == zero

# Counterterm structures follow the local kinetic, mass and vertex operators.
# M=mUV²; deltaZAm multiplies the additive auxiliary mass operator M*A²/2.
# Solve for coefficients of Z=1+a4*deltaZ. No counterterm vertices are generated.
photon_g = uv_poles["photon"].coefficient(gmunu)
photon_pp = uv_poles["photon"].coefficient(ppmunu)
# Unknown order: deltaZpsi, deltaZm, deltaZA, deltaZxi, deltaZe, deltaZAm.
ct_matrix = Matrix.from_linear(
    6,
    6,
    [
        1,
        0,
        0,
        0,
        0,
        0,
        -1,
        -1,
        0,
        0,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        0,
        0,
        0,
        xi - 1,
        1,
        0,
        0,
        1,
        0,
        E("1/2"),
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
    ],
)
rhs = Matrix.vec(
    [
        -uv_poles["electron_p"],
        -uv_poles["electron_m"],
        -photon_g.coefficient(s),
        -xi * photon_pp,
        -uv_poles["vertex"],
        -photon_g.replace(s, zero) / M,
    ]
)
deltas = ct_matrix.solve(rhs)
for row, reference in enumerate(
    [
        -xi / epsilon,
        -3 / epsilon,
        -4 * Nf / (3 * epsilon),
        -4 * Nf / (3 * epsilon),
        2 * Nf / (3 * epsilon),
        zero,
    ]
):
    assert (deltas[row, 0].to_expression() - reference).together() == zero
assert (deltas[0, 0].to_expression() + uv_poles["vertex"]).together() == zero
assert (
    deltas[4, 0].to_expression() + deltas[2, 0].to_expression() / 2
).together() == zero
residual = ct_matrix * deltas - rhs
assert all(residual[row, 0].to_expression().together() == zero for row in range(6))
print(
    "Generated QED one-loop: symbolic gauge, four IBP targets, six counterterms and Ward identities passed",
    solution.stats,
)

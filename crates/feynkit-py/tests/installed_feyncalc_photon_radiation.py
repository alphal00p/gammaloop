"""Generated virtual-photon Born currents and real QCD radiation.

References:
https://feyncalc.github.io/FeynCalcExamples/QED/Tree/Ga-MuAmu
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/Ga-QQbar
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/Ga-QQbarGl
"""

import numpy as np
from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
P = S("gammalooprs::P")
mass, ee, gs = S("UFO::MB", "UFO::ee", "UFO::G")
A, B, C = S("radiation::A", "radiation::B", "radiation::C")
Nc, dA, cof, coad = S("radiation::Nc", "radiation::dA", "spenso::cof", "spenso::coad")
QQ = S("radiation::QQ", is_positive=True)
x1, x2, x3 = S("radiation::x1", "radiation::x2", "radiation::x3")


def calculate(outgoing, kin):
    allowed = [
        v
        for v in model.vertex_rules
        if sorted(v.particles)
        in [
            sorted(["b", "b~", "a"]),
            sorted(["b", "b~", "g"]),
            sorted(["mu-", "mu+", "a"]),
        ]
    ]
    generated = model.generate_diagrams(
        ["a"],
        outgoing,
        vertex_allow=allowed,
        max_vertices=len(outgoing) - 1,
        maximum_bridges=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == len(outgoing) - 1
    ports = S("radiation::i0", "radiation::i1", "radiation::i2", "radiation::i3")
    a, b, c, inv, index, left, right = S(
        "a_", "b_", "c_", "inv_", "index_", "left_", "right_"
    )
    conj, wrap, metric = S("spenso::conj", "radiation::adjoint", "spenso::g")
    terms = []
    for diagram in generated.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        ).replace(S("UFO::MM"), mass)
        for half in diagram.half_edges:
            edge = half.edge.data
            if edge.is_external:
                numerator = numerator.replace(
                    S("gammalooprs::hedge")(half.data, 1), ports[edge.external_index]
                )
        denominator = kin.apply(
            diagram.denominator_expression(dimension=4, in_lmb=True)
            .to_expression()
            .replace(S("gammalooprs::denom")(a, b, c, inv), inv)
        ).expand()
        assert denominator in (E("1"), B, C), denominator
        terms.append(
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
            / denominator
        )
    operator = TensorExpression(sum(terms, E("0")).expand())
    adjoint = (
        operator.dirac_adjoint()
        .expand()
        .simplify_gamma0()
        .to_expression()
        .replace(conj(P(a, b)), P(a, b))
    )
    for real in (mass, ee, gs, A, B, C):
        adjoint = adjoint.replace(conj(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(wrap)
    color_projector = E("1")
    for position, name in enumerate(["a", *outgoing]):
        particle = model.particle(name)
        if particle.color == 1:
            continue
        matches = []
        for slot in operator.interface:
            original = slot.to_expression()
            if original.replace(ports[position], E("0")) == original:
                continue
            closure = metric(
                slot.dual().to_expression(),
                original.replace(ports[position], wrap(ports[position])),
            )
            match = next(
                closure.match(particle.color_sum(left, right), max_level=0), None
            )
            if match is not None:
                matches.append(dict(match))
        assert len(matches) == 1
        color_projector *= particle.color_sum(matches[0][left], matches[0][right])
    generic = (
        (operator.to_expression() * adjoint * color_projector)
        .replace(cof(3, index), cof(Nc, index))
        .replace(coad(8, index), coad(dA, index))
    )
    colored = (
        TensorExpression(generic, cook_indices=CookSettings.indices())
        .simplify_color()
        .to_expression()
        .replace(dA, Nc**2 - 1)
    )
    colored = TensorExpression(colored).to_cof_dimension_invariants()
    print("COLOR PASS", flush=True)
    spins = model.particle(outgoing[0]).spin_sum(
        P(1), wrap(ports[2]), ports[1]
    ) * model.particle(outgoing[1]).spin_sum(P(2), ports[2], wrap(ports[1]))
    traced = (
        TensorExpression(
            colored.to_expression() * spins.replace(S("UFO::MM"), mass),
            cook_indices=CookSettings.indices(),
        )
        .expand()
        .simplify_gamma()
        .expand()
        .simplify_gamma()
    )
    reduced = kin.apply(traced.expand().simplify_metrics().to_dots())
    pieces = [(tensor, coeff.together()) for tensor, coeff in reduced.expand_mink()]
    print("STRUCTURES", len(pieces), flush=True)
    results = {}
    modes = (
        [
            "covariant",
            "quark reference",
            "antiquark reference",
            "gluon Ward",
            "photon Ward",
        ]
        if len(outgoing) == 3
        else ["covariant", "photon Ward"]
    )
    for mode in modes:
        photon = model.particle("a").spin_sum(
            P(0), ports[0], wrap(ports[0]), covariant=True
        )
        gluon = (
            model.particle("g").spin_sum(
                P(3),
                ports[3],
                wrap(ports[3]),
                covariant=mode in ("covariant", "photon Ward", "gluon Ward"),
                reference=P(1)
                if mode == "quark reference"
                else P(2)
                if mode == "antiquark reference"
                else None,
            )
            if len(outgoing) == 3
            else E("1")
        )
        if mode == "gluon Ward":
            gluon = P(3, S("spenso::mink")(4, ports[3])) * P(
                3, S("spenso::mink")(4, wrap(ports[3]))
            )
        if mode == "photon Ward":
            photon = P(0, S("spenso::mink")(4, ports[0])) * P(
                0, S("spenso::mink")(4, wrap(ports[0]))
            )
        contractions = []
        for tensor, coeff in pieces:
            contracted = (
                TensorExpression(
                    (tensor * photon * gluon).expand(),
                    cook_indices=CookSettings.indices(),
                )
                .simplify_metrics()
                .to_dots()
            )
            assert contracted.is_scalar
            contractions.append(coeff * kin.apply(contracted.to_expression()))
        results[mode] = sum(contractions, E("0")).together()
        print("PASS", outgoing, mode, flush=True)
    for mode in [name for name in modes if name.endswith("Ward")]:
        assert results[mode] == E("0")
    for mode in [name for name in modes if name.endswith("reference")]:
        assert (results[mode] - results["covariant"]).together() == E("0")
    print("ALL GAUGES PASS", flush=True)

    return generated, results


kin = hep.Kinematics()
for i, j, value in [
    (1, 1, mass**2),
    (2, 2, mass**2),
    (3, 3, E("0")),
    (1, 2, A / 2),
    (1, 3, B / 2),
    (2, 3, C / 2),
    (0, 0, A + B + C + 2 * mass**2),
    (0, 1, mass**2 + (A + B) / 2),
    (0, 2, mass**2 + (A + C) / 2),
    (0, 3, (B + C) / 2),
]:
    kin = kin.with_scalar_product(P(i), P(j), value)
generated, results = calculate(["b", "b~", "g"], kin)
squared = (
    results["covariant"]
    .replace(A, QQ * (1 - x3))
    .replace(B, QQ * (1 - x2))
    .replace(C, QQ * (1 - x1))
    .together()
)
expected = (
    8
    * ee**2
    * (Nc**2 - 1)
    / 2
    / 9
    * gs**2
    * (
        2
        * QQ
        * mass**2
        * (
            x1**3
            + x1**2 * (x2 + x3 - 5)
            + x1 * (x2**2 - 4 * x2 * x3 + 2 * x3 + 4)
            + x2**3
            + x2**2 * (x3 - 5)
            + 2 * x2 * (x3 + 2)
            - 2 * (x3 + 1)
        )
        - 8 * mass**4 * (x1**2 - 2 * x1 + x2**2 - 2 * x2 + 2)
        + QQ**2
        * (x1 - 1)
        * (x2 - 1)
        * (x1**2 + 2 * x1 * (x3 - 2) + x2**2 + 2 * x2 * (x3 - 2) + 2 * (x3 - 2) ** 2)
    )
    / (QQ**2 * (x1 - 1) ** 2 * (x2 - 1) ** 2)
)
assert (squared - expected).together() == E("0")
massless = squared.replace(mass, 0).replace(x3, 2 - x1 - x2).together()
assert (
    massless
    - 4 * ee**2 * gs**2 * (Nc**2 - 1) / 9 * (x1**2 + x2**2) / ((1 - x1) * (1 - x2))
).together() == E("0")
print("MASSIVE AND MASSLESS REFERENCES PASS", flush=True)

born_kin = hep.Kinematics()
for i, j, value in [
    (0, 0, QQ),
    (1, 1, mass**2),
    (2, 2, mass**2),
    (1, 2, (QQ - 2 * mass**2) / 2),
    (0, 1, QQ / 2),
    (0, 2, QQ / 2),
]:
    born_kin = born_kin.with_scalar_product(P(i), P(j), value)
born_results = {}
for names, color_charge in ((["b", "b~"], Nc / 9), (["mu-", "mu+"], E("1"))):
    born_generated, born_sums = calculate(names, born_kin)
    born_square = born_sums["covariant"]
    assert (
        born_square - 4 * ee**2 * color_charge * (QQ + 2 * mass**2)
    ).together() == E("0")
    # The virtual-photon current is contracted with -g_mu_nu without averaging,
    # matching the gallery normalization. This is not an on-shell photon decay.
    width = (
        born_square
        * 4
        * Symbol.PI
        * born_kin.two_body_phase_space(P(1), P(2))
        / born_kin.flux(P(0))
    )
    massless_width = width.replace(mass, 0).together()
    assert (
        massless_width - ee**2 * color_charge * QQ.sqrt() / (4 * Symbol.PI)
    ).together() == E("0")
    born_results[names[0]] = (born_generated, born_square, width, massless_width)
assert (born_results["b"][3] / born_results["mu-"][3] - Nc / 9).together() == E("0")
print("PASS: generated massive quark/lepton Born currents, massless widths and R ratio")

massless_kin = hep.Kinematics()
for i, j, value in [
    (0, 0, QQ),
    (1, 1, E("0")),
    (2, 2, E("0")),
    (3, 3, E("0")),
    (1, 2, QQ * (x1 + x2 - 1) / 2),
    (1, 3, QQ * (1 - x2) / 2),
    (2, 3, QQ * (1 - x1) / 2),
]:
    massless_kin = massless_kin.with_scalar_product(P(i), P(j), value)
phase = massless_kin.three_body_phase_space(P(1), P(2), P(3)).together()
assert (phase - 1 / (128 * Symbol.PI**3 * QQ)).together() == E("0")
# The transformation from pair invariant masses to energy fractions has Jacobian Q^4.
alpha_s = S("radiation::alpha_s")
distribution = (
    (massless * phase * QQ**2 / massless_kin.flux(P(0)) / born_results["b"][3])
    .replace(gs**2, 4 * Symbol.PI * alpha_s)
    .together()
)
cf = (Nc**2 - 1) / (2 * Nc)
shape = (x1**2 + x2**2) / ((1 - x1) * (1 - x2))
assert (distribution - alpha_s * cf / (2 * Symbol.PI) * shape).together() == E("0")
# The full generated result must reproduce the independently audited soft limit.
lam = S("radiation::lambda")
soft_ratio = results["covariant"].replace(mass, 0) / (4 * ee**2 * Nc * A / 9)
soft_limit = (
    soft_ratio.replace(B, lam**2 * B)
    .replace(C, lam**2 * C)
    .series(lam, 0, -4)
    .to_expression()
    .replace(lam, 1)
)
assert (soft_limit - 4 * gs**2 * cf * A / (B * C)).together() == E("0")
print("PASS: native three-body phase space and Born-normalized Dalitz distribution")

# Positive invariants y1=1-x1 and y2=1-x2 obey y1+y2<=1. The cut is
# y1,y2>=beta, 0<beta<1/2, excluding both soft/collinear poles.
y1, y2, beta = S("radiation::y1", "radiation::y2", "radiation::beta", is_positive=True)
kernel = shape.replace(x1, 1 - y1).replace(x2, 1 - y2).together()
argument = S("radiation::argument_")
inner_primitive = kernel.integrate(y2).replace(
    Symbol.LOG(argument), lambda match: match[argument].expand().log()
)
assert (inner_primitive.derivative(y2) - kernel).together() == E("0")
inner_cut = inner_primitive.replace(y2, 1 - y1) - inner_primitive.replace(y2, beta)
# The triangular domain is symmetric under y1<->y2. Integrating the two
# symmetric numerator terms therefore gives twice the first term's integral.
outer_integrand = 2 * (1 - y1) ** 2 / y1 * ((1 - y1).log() - beta.log())
outer_primitive = (
    outer_integrand.integrate(y1)
    .replace(
        Symbol.POLYLOG(2, argument), lambda match: match[argument].expand().polylog(2)
    )
    .replace(Symbol.LOG(argument), lambda match: match[argument].expand().log())
)
assert (outer_primitive.derivative(y1) - outer_integrand).together() == E("0")
integrated_kernel = outer_primitive.replace(y1, 1 - beta) - outer_primitive.replace(
    y1, beta
)
# A real-branch form follows from Euler's dilogarithm reflection identity.
closed = (
    2 * beta.log() ** 2
    + (3 - 4 * beta + beta**2) * (beta.log() - (1 - beta).log())
    + E("5/2")
    - 5 * beta
    + 4 * beta.polylog(2)
    - Symbol.PI**2 / 3
)
# Expand the analytic coefficients while retaining log(beta) symbolically.
log_beta = S("radiation::log_beta")
leading = (
    closed.replace(beta.log(), log_beta)
    .series(beta, 0, 0)
    .to_expression()
    .replace(log_beta, beta.log())
)
expected_leading = 2 * beta.log() ** 2 + 3 * beta.log() + E("5/2") - Symbol.PI**2 / 3
assert (leading - expected_leading).expand() == E("0")
# Leibniz differentiation of the cut triangle provides an independent exact
# check of the closed form, without relying on reflection simplification.
boundary = inner_cut.replace(y1, beta)
assert (closed.derivative(beta) + 2 * boundary).together() == E("0")
nodes, weights = np.polynomial.legendre.leggauss(96)
for cutoff in (0.01, 0.05, 0.1, 0.25, 0.4):
    first = cutoff + (nodes + 1) * (1 - 2 * cutoff) / 2
    total = 0.0
    for u, wu in zip(first, weights):
        second = cutoff + (nodes + 1) * (1 - u - cutoff) / 2
        values = ((1 - u) ** 2 + (1 - second) ** 2) / (u * second)
        total += wu * (1 - u - cutoff) / 2 * float(np.dot(weights, values))
    total *= (1 - 2 * cutoff) / 2
    exact = complex(closed.evaluate({beta: cutoff}))
    automatic = complex(integrated_kernel.evaluate({beta: cutoff}))
    assert abs(exact.imag) < 1e-12 and abs(automatic.imag) < 1e-12
    assert abs(exact.real - total) < 2e-10
    assert abs(automatic.real - total) < 2e-10
    for nc in (2, 3, 5):
        rate = (alpha_s * cf / (2 * Symbol.PI) * closed).evaluate(
            {beta: cutoff, Nc: nc, alpha_s: 0.118}
        )
        assert complex(rate).real > 0
assert abs(complex(closed.evaluate({beta: 0.5}))) < 1e-12
print(
    "PASS: both automatic primitives, cut-boundary identity, IR logarithms and 15 independent quadratures"
)

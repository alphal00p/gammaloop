"""Generated massive quark scattering, color interference and angular-cut rates.

References: FeynCalc QCD/Tree QiQj-QiQj, QiQjbar-QiQjbar,
QiQi-QiQi and QiQibar-QiQibar. All four massive SU(N) expressions,
interference terms and massless limits are checked before phase-space integration.
"""

import numpy as np
from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
P = S("gammalooprs::P")
s, t, u, m, M, gs = S(
    "qscatter::s", "qscatter::t", "qscatter::u", "UFO::MB", "UFO::MT", "UFO::G"
)
Nc, dA, cof, coad = S("qscatter::Nc", "qscatter::dA", "spenso::cof", "spenso::coad")
a, b, c, inverse, wave, rep, index = S(
    "a_", "b_", "c_", "inverse_", "wave_", "rep_", "index_"
)
ports = S("qscatter::i0", "qscatter::i1", "qscatter::i2", "qscatter::i3")
conj, wrapped = S("spenso::conj", "qscatter::adjoint")
flavors = [model.particle_by_pdg(i) for i in (5, 6)]
gluon = model.particle_by_pdg(21)
vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles)
    in [sorted([q.name, q.antiname, gluon.name]) for q in flavors]
]
assert len(vertices) == 2
results = {}
for name, pdgs in [
    ("qq_prime", [5, 6, 5, 6]),
    ("qaq_prime", [5, -6, 5, -6]),
    ("qq", [5, 5, 5, 5]),
    ("qaq", [5, -5, 5, -5]),
]:
    same = name in ("qq", "qaq")
    masses = [m**2, (m if same else M) ** 2] * 2
    kin = hep.Kinematics.mandelstam([P(i) for i in range(4)], masses, [s, t, u])
    generation = model.process(
        pdgs[:2], pdgs[2:], vertex_allow=vertices
    ).generate_diagrams(
        max_vertices=2, maximum_bridges=None, numerator_grouping=None, progress=None
    )
    assert len(generation.diagrams) == (2 if same else 1)
    operators, denominators = [], []
    for diagram in generation.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
        for edge in diagram.external_edges:
            match = dict(
                next(
                    diagram.projector_expression().match(
                        wave(edge.id, rep(4, index)), max_level=0
                    )
                )
            )
            numerator = numerator.replace(match[index], ports[edge.external_index])
        denominator = kin.apply(
            diagram.denominator_expression(dimension=4, in_lmb=True)
            .to_expression()
            .replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
        ).expand()
        denominators.append(denominator)
        operators.append(
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
            / denominator
        )
    assert set(denominators) == (
        {s, t} if name == "qaq" else {t, u} if name == "qq" else {t}
    )
    amplitude = sum(operators, E("0"))
    operator = TensorExpression(amplitude.expand())
    assert len(operator.structure.slots) == 8
    adjoints = []
    for term in operators:
        adjoint = (
            TensorExpression(term).dirac_adjoint(preserve_indices=True).to_expression()
        )
        for real in (s, t, u, m, M, gs):
            adjoint = adjoint.replace(conj(real), real)
        adjoints.append(TensorExpression(adjoint).wrap_indices(wrapped))
    colors = E("1")
    spins = E("1")
    initial_colors = 1
    for position, pdg in enumerate(pdgs):
        particle = model.particle_by_pdg(pdg)
        if position < 2:
            initial_colors *= abs(particle.color)
        # Signed particle representations choose quark/antiquark duals.
        # Incoming color ports close in the reverse order from outgoing ports.
        ci, cj = (
            (wrapped(ports[position]), ports[position])
            if position < 2
            else (ports[position], wrapped(ports[position]))
        )
        colors *= particle.color_sum(ci, cj, average=position < 2)
        column = (position < 2) == (pdg > 0)
        i, j = (
            (ports[position], wrapped(ports[position]))
            if column
            else (wrapped(ports[position]), ports[position])
        )
        spins *= particle.spin_sum(P(position), i, j, average=position < 2)
    evaluated = []
    for raw in (
        amplitude * sum(adjoints, E("0")),
        sum((x * y for x, y in zip(operators, adjoints, strict=True)), E("0")),
    ):
        # Promote the model's initial SU(3) color average to SU(N). The
        # adjoint dimension remains a symbol until color contraction finishes.
        generic = (
            (raw * colors * initial_colors / Nc**2)
            .replace(cof(3, index), cof(Nc, index))
            .replace(coad(8, index), coad(dA, index))
        )
        colored = (
            TensorExpression(generic, cook_indices=CookSettings.indices())
            .simplify_color()
            .to_expression()
            .replace(dA, Nc**2 - 1)
        )
        colored = (
            TensorExpression(colored).to_cof_dimension_invariants().to_expression()
        )
        scalar = (
            TensorExpression(colored * spins, cook_indices=CookSettings.indices())
            .expand()
            .simplify_gamma()
            .expand()
            .simplify_gamma()
            .simplify_metrics()
            .to_dots()
        )
        assert scalar.is_scalar
        value = (
            kin.apply(scalar.to_expression()).replace(u, sum(masses) - s - t).together()
        )
        evaluated.append(value)
    squared, diagonal = evaluated
    if same:
        x, y, z = (t, u, s) if name == "qq" else (s, t, u)
        expected = (
            (Nc**2 - 1)
            * gs**4
            * (
                -4 * m**2 * (Nc * (x**3 + y**3) - 2 * s * t * u)
                + 4 * m**4 * (Nc * (x**2 + y**2) - 3 * x * y)
                + Nc * (x**4 + x**3 * y + x**2 * y**2 + x * y**3 + y**4)
                - z**2 * x * y
            )
            / (Nc**3 * x**2 * y**2)
        )
        interference = (
            -(Nc**2 - 1) * gs**4 * (z**2 - 8 * m**2 * z + 12 * m**4) / (Nc**3 * x * y)
        )
        massless_expected = (
            (Nc**2 - 1)
            * gs**4
            / (2 * Nc**2)
            * (
                (s**2 + u**2) / t**2
                + ((s**2 + t**2) / u**2 if name == "qq" else (t**2 + u**2) / s**2)
            )
        )
        massless_expected -= (Nc**2 - 1) * gs**4 * z**2 / (Nc**3 * x * y)
    else:
        expected = (
            (Nc**2 - 1)
            * gs**4
            * (
                -4 * M**2 * (u - m**2)
                + 2 * M**4
                + 2 * m**4
                - 4 * u * m**2
                + t**2
                + 2 * t * u
                + 2 * u**2
            )
            / (2 * Nc**2 * t**2)
        )
        interference = E("0")
        massless_expected = (Nc**2 - 1) * gs**4 * (s**2 + u**2) / (2 * Nc**2 * t**2)
    residual = (squared - expected.replace(u, sum(masses) - s - t)).together()
    assert residual == 0
    assert (
        squared - diagonal - interference.replace(u, sum(masses) - s - t)
    ).together() == 0
    massless = squared.replace(m, E("0")).replace(M, E("0"))
    assert (massless - massless_expected.replace(u, -s - t)).together() == 0
    if name == "qq":
        assert (squared - squared.replace(t, 4 * m**2 - s - t)).together() == 0
    phase = kin.two_body_phase_space(P(2), P(3)) / kin.flux(P(0), P(1))
    assert (phase - 1 / (64 * Symbol.PI**2 * s)).together() == 0
    results[name] = {
        "squared": squared,
        "diagonal": diagonal,
        "interference": squared - diagonal,
        "diagrams": generation.diagrams,
        "kinematics": kin,
        "same_flavor": same,
    }
    print(
        name,
        "massive SU(N), interference, massless and native phase-space checks passed",
        flush=True,
    )

# Carry the generated massive result through the native elastic flux and phase
# space. Full angular integration needs a cut because of massless gluon exchange.
z, cutoff, alpha_s, C = S(
    "qscatter::z", "qscatter::cutoff", "qscatter::alpha_s", "qscatter::C"
)
template = sum((C(i) * z**i for i in range(7)), E("0")) / (1 - z**2) ** 2
primitive_template = template.integrate(z)
assert (primitive_template.derivative(z) - template).together() == 0
rate_checks = []
nodes, weights = np.polynomial.legendre.leggauss(128)
for name, result in results.items():
    same = result["same_flavor"]
    second_mass = m if same else M
    kallen = (s - m**2 - second_mass**2) ** 2 - 4 * m**2 * second_mass**2
    angle_t = -kallen * (1 - z) / (2 * s)
    symmetry = E("1/2") if name == "qq" else E("1")
    prefactor = symmetry / (32 * Symbol.PI * s)
    density = (
        result["squared"]
        .replace(t, angle_t)
        .replace(gs**4, (4 * Symbol.PI * alpha_s) ** 2)
        * prefactor
    ).together()
    coefficient_polynomial = (density * (1 - z**2) ** 2).together().expand()
    coefficients = {
        monomial.to_polynomial(vars=[z]).degree(z): coefficient
        for monomial, coefficient in coefficient_polynomial.coefficient_list(z)
    }
    assert set(coefficients) <= set(range(7))
    primitive = primitive_template.replace_multiple(
        [Replacement(C(i), coefficients.get(i, E("0"))) for i in range(7)]
    )
    assert (primitive.derivative(z) - density).together() == 0
    # On -1<z<1, replacing log(z-1) by log(1-z) changes only an irrelevant
    # additive imaginary constant. This gives a manifestly real primitive.
    primitive = primitive.replace((z - 1).log(), (1 - z).log())
    cut_rate = primitive.replace(z, cutoff) - primitive.replace(z, -cutoff)
    assert cut_rate.replace(cutoff, 0).together() == 0
    if name == "qq":
        # One forward quark per event equals half the full labeled phase space.
        residual = (
            primitive.replace(z, cutoff) - primitive.replace(z, 0) - cut_rate / 2
        ).together()
        # The real-interval identity can retain atanh(-cutoff) symbolically.
        # Equal derivatives and the value at zero establish the exact equality.
        assert residual.derivative(cutoff).together() == 0
        assert residual.replace(cutoff, 0).together() == 0
    result.update(
        density=density, cut_rate=cut_rate, angle_t=angle_t, symmetry=symmetry
    )
    for nv in (2, 3, 5):
        for fraction, other_fraction in [(0.0, 0.0), (0.1, 0.15), (0.2, 0.25)]:
            for sv in (4.0, 25.0):
                for cv in (0.2, 0.5, 0.8):
                    values = {
                        Nc: nv,
                        m: fraction * np.sqrt(sv),
                        M: other_fraction * np.sqrt(sv),
                        s: sv,
                        alpha_s: 0.118,
                    }
                    integrated = complex(cut_rate.evaluate({**values, cutoff: cv}))
                    quadrature = cv * sum(
                        weight * complex(density.evaluate({**values, z: cv * node}))
                        for node, weight in zip(nodes, weights, strict=True)
                    )
                    assert abs(integrated.imag) < 1e-11
                    assert integrated.real > 0
                    error = abs(integrated - quadrature)
                    assert error < 2e-11 * max(1, abs(integrated)), (
                        name,
                        nv,
                        fraction,
                        sv,
                        cv,
                        error,
                    )
                    rate_checks.append((name, nv, fraction, sv, cv, error))
    print(name, "massive angular-cut rates and event counting passed", flush=True)
print(len(rate_checks), "independent phase-space quadratures passed", flush=True)

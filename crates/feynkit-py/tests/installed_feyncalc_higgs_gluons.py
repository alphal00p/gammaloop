"""Generated top-quark H -> gg triangle, native tensor/IBP and OneLOop.

Reference: https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/H-GlGl
The gallery uses CreateFeynAmp PreFactor -> -1. Native UFO amplitudes retain
all graph/Feynman-rule phases, so their overall sign is opposite to that page.
"""

from math import asin, asinh, log, pi, sqrt
from pathlib import Path

from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import ColorSimplifySettings, TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
D, s, m, gs, y, mu, nu, ca, cb, eps = S(
    "hgg::D",
    "hgg::s",
    "UFO::MT",
    "UFO::G",
    "hgg::y",
    "hgg::mu",
    "hgg::nu",
    "hgg::ca",
    "hgg::cb",
    "hgg::eps",
)
K, P, mink, cof, coad, metric = S(
    "gammalooprs::K",
    "gammalooprs::P",
    "spenso::mink",
    "spenso::cof",
    "spenso::coad",
    "spenso::g",
)
index, wave, arguments = S("index_", "wave_", "arguments__")
Nc, dA = S("hgg::Nc", "hgg::dA")
# The SM UFO stores yt=sqrt(2)*m/v separately from the propagator pole mass.
# Define y=yt/sqrt(2); the comparison later imposes y=m/v explicitly.
top, higgs, gluon = (model.particle_by_pdg(pdg) for pdg in (6, 25, 21))
yukawa_vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles) == sorted([top.antiname, top.name, higgs.name])
]
gluon_vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles) == sorted([top.antiname, top.name, gluon.name])
]
assert len(yukawa_vertices) == len(gluon_vertices) == 1
yukawa_tree_result = model.generate_diagrams(
    [higgs],
    [top, top.antiparticle],
    max_vertices=1,
    vertex_allow=yukawa_vertices,
    numerator_grouping=None,
    progress=None,
)
assert len(yukawa_tree_result.diagrams) == 1
yukawa_tree = yukawa_tree_result.diagrams[0]
yukawa_tree_kernel = model.expand_couplings(
    yukawa_tree.numerator_expression().to_expression()
).replace(S("UFO::yt") * E("1/2").sqrt(), y)
assert yukawa_tree.overall_factor_expression(evaluate=True) == E("1")
assert yukawa_tree.numerator_prefactor_expression() == E("1")
# Strip only the generated spin/color identity to establish the tree phase.
identities = S("left_", "right_")
yukawa_tree_coupling = yukawa_tree_kernel.replace(metric(*identities), E("1"))
assert (yukawa_tree_coupling + Symbol.I * y).expand() == E("0")

result = model.generate_diagrams(
    [higgs],
    [gluon, gluon],
    loops=1,
    max_vertices=3,
    maximum_bridges=0,
    vertex_allow=gluon_vertices + yukawa_vertices,
    numerator_grouping=None,
    progress=None,
)
assert len(result.diagrams) == 2
# P0=k1+k2 is incoming Higgs momentum and P1=k1 the first outgoing gluon.
kinematics = (
    hep.Kinematics(D, momenta=[K(0), P(0), P(1)])
    .with_scalar_product(P(0), P(0), s)
    .with_scalar_product(P(1), P(1), E("0"))
    .with_scalar_product(P(0), P(1), s / 2)
)
reducer = (
    hep.TensorReducer(D)
    .with_integrated_vector(K(0, mink(D)))
    .with_external_vector(P(0, mink(D)))
    .with_external_vector(P(1, mink(D)))
)
family = result.diagrams[0].integral_family(kinematics=kinematics)
coordinates = S("hgg::d0", "hgg::d1", "hgg::d2")
Qg, Q11, Q12, Q21, Q22 = basis = S(
    "hgg::Qg", "hgg::Q11", "hgg::Q12", "hgg::Q21", "hgg::Q22"
)
p_mu, p_nu, q_mu, q_nu = S("hgg::pmu", "hgg::pnu", "hgg::qmu", "hgg::qnu")
color = gluon.color_sum(ca, cb).replace(coad(8, index), coad(dA, index))
terms_by_diagram, targets, routed_polynomials = [], set(), []
for orientation, diagram in enumerate(result.diagrams):
    assert diagram.overall_factor_expression(evaluate=True) == E("-1")
    assert diagram.numerator_prefactor_expression() == E("1")
    numerator = model.expand_couplings(
        diagram.numerator_expression(in_lmb=True).to_expression()
    ).replace(S("UFO::yt") * E("1/2").sqrt(), y)
    for edge in diagram.external_edges:
        if edge.external_index == 0:
            continue
        port = dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, mink(4, index)), max_level=0
                )
            )
        )[index]
        numerator = numerator.replace(
            mink(4, port), mink(D, (mu, nu)[edge.external_index - 1])
        )
        numerator = numerator.replace(
            coad(8, port), coad(dA, (ca, cb)[edge.external_index - 1])
        )
    numerator = (
        numerator.replace(mink(4, index), mink(D, index))
        .replace(cof(3, index), cof(Nc, index))
        .replace(coad(8, index), coad(dA, index))
    )
    # Shared Idenso evaluates Tr(Ta Tb)=1/2 delta_ab, before summing the two
    # orientations. Promote Lorentz dimension before Dirac trace; tr(1)=4.
    trace = (
        TensorExpression(numerator.expand())
        .simplify_gamma()
        .expand()
        .simplify_color(ColorSimplifySettings(substitute_cof_dimension_invariants=True))
        .simplify_metrics()
        .to_dots()
        .to_expression()
    )
    trace = (
        trace
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    ).expand()
    # The complete trace vanishes for a massless quark even at fixed Yukawa
    # coupling; check chirality before dividing out the explicit mass factor.
    assert trace.replace(m, E("0")).expand() == E("0")
    assert (trace.coefficient(color) * color - trace).expand() == E("0")
    trace = trace.coefficient(color) / (gs**2 * y * m)
    mapping = diagram.integral_family(kinematics=kinematics).mapping_to(
        family, [(-1 if orientation else 1) * K(0)]
    )
    assert mapping is not None
    # IntegralMapping acts on scalar numerators: reduce free loop indices first.
    reduced = mapping.apply(kinematics.apply(reducer.reduce(trace))).expand()
    reduced = reduced.replace_multiple(
        [
            Replacement(metric(mink(D, mu), mink(D, nu)), Qg),
            Replacement(P(0, mink(D, mu)), p_mu + q_mu),
            Replacement(P(0, mink(D, nu)), p_nu + q_nu),
            Replacement(P(1, mink(D, mu)), p_mu),
            Replacement(P(1, mink(D, nu)), p_nu),
        ]
    ).expand()
    reduced = reduced.replace_multiple(
        [
            Replacement(p_mu * p_nu, Q11),
            Replacement(p_mu * q_nu, Q12),
            Replacement(q_mu * p_nu, Q21),
            Replacement(q_mu * q_nu, Q22),
        ]
    )
    polynomial = family.rewrite_numerator(reduced, coordinates).expand()
    routed_polynomials.append(polynomial)
    terms = []
    for monomial, coefficient in polynomial.coefficient_list(*coordinates):
        powers = tuple(
            1 - monomial.to_polynomial(vars=coordinates).degree(label)
            for label in coordinates
        )
        assert not coefficient.matches(K(arguments))
        assert not coefficient.matches(P(arguments))
        assert all(
            coefficient.derivative(label).expand() == E("0") for label in coordinates
        )
        targets.add(powers)
        terms.append((powers, coefficient))
    terms_by_diagram.append(terms)
assert (routed_polynomials[0] - routed_polynomials[1]).together() == E("0")
solution = hep.IBPFamily(family, name="higgs_gluon_triangle").reduce_laporta(
    [list(t) for t in sorted(targets)], max_depth=2
)
assert {tuple(p) for p in solution.residuals} == {
    (0, 0, 1),
    (0, 1, 0),
    (0, 1, 1),
    (1, 0, 0),
    (1, 1, 1),
}
integral, a0, b0, c0 = S("hgg::I", "hgg::A0", "hgg::B0", "hgg::C0")
# One-line pinches are shifted equal-mass tadpoles; the two-line pinch has
# invariant s, and the three-line master has external invariants (0,0,s).
for powers in ((0, 0, 1), (0, 1, 0), (1, 0, 0)):
    assert family.sector(powers).find_mapping(family.sector([0, 1, 0])) is not None
# Symanzik polynomials independently identify the mass/invariant arguments
# of the scalar B0(s;m^2,m^2) and C0(0,0,s;m^2,m^2,m^2) masters.
x0, x1, x2 = S("hgg::x0", "hgg::x1", "hgg::x2")
triangle_U, triangle_F = family.symanzik([x0, x1, x2])
assert (triangle_U - x0 - x1 - x2).expand() == E("0")
assert (triangle_F - m**2 * triangle_U**2 + s * x1 * x2).expand() == E("0")
bubble_U, bubble_F = family.sector([0, 1, 1]).symanzik([x1, x2])
assert (bubble_U - x1 - x2).expand() == E("0")
assert (bubble_F - m**2 * bubble_U**2 + s * x1 * x2).expand() == E("0")
masters = [Replacement(integral(*p), a0) for p in ((0, 0, 1), (0, 1, 0), (1, 0, 0))]
masters += [Replacement(integral(0, 1, 1), b0), Replacement(integral(1, 1, 1), c0)]
integrated_by_diagram = []
for terms in terms_by_diagram:
    reduction = sum(
        (
            coefficient * solution.reduce(list(powers), integral=integral)
            for powers, coefficient in terms
        ),
        E("0"),
    )
    integrated_by_diagram.append(reduction.replace_multiple(masters).together())
assert (integrated_by_diagram[0] - integrated_by_diagram[1]).together() == E("0")
integrated = sum(integrated_by_diagram, E("0")).expand()
# Truncating scalar masters at finite order is safe: their rational
# coefficients have no hidden (D-4) poles.
for master in (a0, b0, c0):
    assert (
        integrated.coefficient(master)
        .together()
        .replace(D, 4 - 2 * eps)
        .series(eps, 0, -1)
        .to_expression()
    ) == E("0")
assert integrated.coefficient(Q11).together() == E("0")
assert integrated.coefficient(Q22).together() == E("0")
# On-shell open tensors can retain k1_mu*k2_nu. It is transverse to both
# massless momenta and vanishes against their own physical polarizations.
physical = integrated.replace(Q12, E("0")).together()
transverse = s * Qg / 2 - Q21
form_factor = (physical.expand().coefficient(Qg) * 2 / s).together()
assert (physical - form_factor * transverse).together() == E("0")
expected_D = -4 / s * (2 * (4 - D) / (D - 2) * b0 + (8 * m**2 / (D - 2) - s) * c0)
assert (form_factor - expected_D).together() == E("0")
# Keep D dependence until after supplying the scalar pole residues. The
# finite rational 2 is produced by (D-4)*B0; setting D=4 first loses it.
Af, Bf, Cf = S("hgg::Af", "hgg::Bf", "hgg::Cf")
laurent = (
    integrated.replace(D, 4 - 2 * eps)
    .replace_multiple(
        [
            Replacement(a0, m**2 / eps + Af),
            Replacement(b0, 1 / eps + Bf),
            Replacement(c0, Cf),
        ]
    )
    .series(eps, 0, 0)
    .to_expression()
    .expand()
)
assert laurent.coefficient(eps**-1).together() == E("0")
assert laurent.coefficient(eps**-2).together() == E("0")
finite = dict(laurent.coefficient_list(eps))[E("1")].together().expand()
finite_factor = (finite.coefficient(Qg) * 2 / s).together()
assert (finite_factor + 4 / s * (2 + (4 * m**2 - s) * Cf)).together() == E("0")
# Scalar masters restore i/(16*pi^2). Thus native iM has coefficient
# -i*gs^2*y*m/(4*pi^2*s) [2+(4m^2-s)C0] multiplying delta_ab*T_mu_nu.
normalized = (-3 * m**2 * finite_factor / 4).together()
assert (normalized - 3 * m**2 / s * (2 + (4 * m**2 - s) * Cf)).together() == E("0")
small_s_triangle = (
    -1 / (2 * m**2) - s / (24 * m**4) - s**2 / (180 * m**6) - s**3 / (1120 * m**8)
)
heavy_series = normalized.replace(Cf, small_s_triangle).series(s, 0, 2).to_expression()
assert (heavy_series - 1 - 7 * s / (120 * m**2) - s**2 / (168 * m**4)).expand() == E(
    "0"
)

# Shared spin/color contractions establish the unaveraged tensor norm.
tensor = s * metric(mink(D, mu), mink(D, nu)) / 2 - (
    P(0, mink(D, mu)) - P(1, mink(D, mu))
) * P(1, mink(D, nu))
exchanged = tensor.replace_multiple(
    [
        Replacement(P(0, mink(D, mu)), P(0, mink(D, nu))),
        Replacement(P(1, mink(D, mu)), P(0, mink(D, nu)) - P(1, mink(D, nu))),
        Replacement(P(1, mink(D, nu)), P(0, mink(D, mu)) - P(1, mink(D, mu))),
    ]
)
assert (tensor - exchanged).expand() == E("0")
for contraction in (P(1, mink(D, mu)), P(0, mink(D, nu)) - P(1, mink(D, nu))):
    ward = (
        TensorExpression((tensor * contraction).expand())
        .simplify_metrics()
        .to_dots()
        .to_expression()
    )
    assert kinematics.apply(ward).expand() == E("0")
norm = kinematics.apply(
    TensorExpression((tensor**2).expand()).simplify_metrics().to_dots().to_expression()
)
assert (norm - (D - 2) * s**2 / 4).expand() == E("0")
color_norm = TensorExpression(color**2).simplify_metrics().to_expression()
assert color_norm == dA
# Also evaluate the physical axial polarization sums, taking each gluon as
# the other's null reference. Ward identities imply agreement with -g sums.
rho, sigma = S("hgg::rho", "hgg::sigma")
physical_kinematics = (
    hep.Kinematics(momenta=[P(1), P(2)])
    .with_scalar_product(P(1), P(1), E("0"))
    .with_scalar_product(P(2), P(2), E("0"))
    .with_scalar_product(P(1), P(2), s / 2)
)
physical_tensor = s * metric(mink(4, mu), mink(4, nu)) / 2 - P(2, mink(4, mu)) * P(
    1, mink(4, nu)
)
conjugate_tensor = physical_tensor.replace_multiple(
    [Replacement(mu, rho), Replacement(nu, sigma)]
)
polarization_sum = gluon.spin_sum(P(1), mu, rho, reference=P(2)) * gluon.spin_sum(
    P(2), nu, sigma, reference=P(1)
)
physical_norm = physical_kinematics.apply(
    TensorExpression((physical_tensor * conjugate_tensor * polarization_sum).expand())
    .simplify_metrics()
    .to_dots()
    .to_expression()
)
assert (physical_norm - s**2 / 2).together() == E("0")
# Contract the actual full finite open tensor, including its Q12 coefficient,
# against both physical polarization projectors. Their difference is exactly
# zero without setting the coefficient of Q12 to zero by hand.
full_finite_tensor = finite.replace_multiple(
    [
        Replacement(Qg, metric(mink(4, mu), mink(4, nu))),
        Replacement(Q12, P(1, mink(4, mu)) * P(2, mink(4, nu))),
        Replacement(Q21, P(2, mink(4, mu)) * P(1, mink(4, nu))),
    ]
)
physical_projection_residual = physical_kinematics.apply(
    TensorExpression(
        (
            (full_finite_tensor - finite_factor * physical_tensor) * polarization_sum
        ).expand()
    )
    .simplify_metrics()
    .to_dots()
    .to_expression()
).together()
assert physical_projection_residual == E("0")
# y=m/v and gs^2=4*pi*alpha_s. Color dimension is specialized after contraction.
alpha_s, vev, A, Abar = S("hgg::alpha_s", "hgg::vev", "hgg::A", "hgg::Abar")
MH = S("hgg::MH", is_positive=True)
squared = (
    physical_norm * color_norm * (alpha_s / (3 * Symbol.PI * vev)) ** 2 * A * Abar
).replace(dA, Nc**2 - 1)
decay_kinematics = (
    hep.Kinematics()
    .with_scalar_product(P(0), P(0), MH**2)
    .with_scalar_product(P(1), P(1), E("0"))
    .with_scalar_product(P(2), P(2), E("0"))
    .with_scalar_product(P(1), P(2), MH**2 / 2)
)
phase_space = decay_kinematics.two_body_phase_space(P(1), P(2)).expand()
flux = decay_kinematics.flux(P(0))
assert (phase_space - 1 / (32 * Symbol.PI**2)).together() == E("0")
assert flux == 2 * MH
# Integrate dOmega=4*pi and include 1/2! for identical final gluons.
width = (squared.replace(s, MH**2) * 4 * Symbol.PI * phase_space / flux / 2).together()
assert (
    width - (Nc**2 - 1) * alpha_s**2 * MH**3 * A * Abar / (576 * Symbol.PI**3 * vev**2)
).together() == E("0")

# Shared OneLOop continuation is checked against the independent analytic
# Feynman-parameter result. Above threshold its -i0 boundary gives Im(C0)<0;
# the physical squared form factor uses the complex modulus, never A^2.
scale2 = S("hgg::scale2")
a_coefficients = oneloop.master_coefficients(S("oneloopmaster::A0")(m**2, scale2))
b_coefficients = oneloop.master_coefficients(
    S("oneloopmaster::B0")(s, m**2, m**2, scale2)
)
c_coefficients = oneloop.master_coefficients(
    S("oneloopmaster::C0")(0, 0, s, m**2, m**2, m**2, scale2)
)
finite_oneloop = normalized.replace(Cf, c_coefficients[0])
numeric_points = [
    (-0.1, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
    (-4.0, 1.0, 1.0),
    (-40.0, 1.0, 3.0),
    (-3.0, 2.0, 7.0),
    (0.001, 1.0, 1.0),
    (0.04, 1.0, 0.25),
    (0.4, 1.0, 1.0),
    (2.0, 1.0, 3.0),
    (3.99, 1.0, 1.0),
    (4.0, 1.0, 1.0),
    (4.01, 1.0, 1.0),
    (8.0, 1.0, 1.0),
    (100.0, 1.0, 7.0),
    (10000.0, 1.0, 1.0),
    (3.0, 2.0, 7.0),
    (16.0, 2.0, 3.0),
    (20.0, 2.0, 7.0),
    (64.0, 0.5, 2.0),
]
numeric_values = []
for invariant, mass, scale in numeric_points:
    point = {s: invariant, m: mass, scale2: scale}
    if invariant < 0:
        c_reference = 2 * asinh(sqrt(-invariant) / (2 * mass)) ** 2 / invariant
    elif invariant <= 4 * mass**2:
        c_reference = -2 * asin(sqrt(invariant) / (2 * mass)) ** 2 / invariant
    else:
        beta = sqrt(1 - 4 * mass**2 / invariant)
        c_reference = (log((1 + beta) / (1 - beta)) - 1j * pi) ** 2 / (2 * invariant)
    c_numeric = complex(c_coefficients[0].evaluate(point))
    assert abs(c_numeric - c_reference) < 2e-12, (point, c_numeric, c_reference)
    assert abs(complex(c_coefficients[1].evaluate(point))) < 1e-12
    assert abs(complex(c_coefficients[2].evaluate(point))) < 1e-12
    assert abs(complex(a_coefficients[1].evaluate(point)) - mass**2) < 1e-12
    assert abs(complex(b_coefficients[1].evaluate(point)) - 1) < 1e-12
    # This finite triangle is independent of the dimensional scale.
    scale_shifted = complex(c_coefficients[0].evaluate({**point, scale2: 7 * scale}))
    assert abs(c_numeric - scale_shifted) < 2e-12
    value = complex(finite_oneloop.evaluate(point))
    reference = 3 * mass**2 / invariant * (2 + (4 * mass**2 - invariant) * c_reference)
    # The heavy limit subtracts O(1) terms before division by s/m^2.
    tolerance = 3e-12 * max(1.0, mass**2 / abs(invariant))
    assert abs(value - reference) < tolerance, (point, value, reference)
    if invariant > 4 * mass**2:
        assert c_numeric.imag < 0 and value.imag > 0
    else:
        assert abs(value.imag) < tolerance
    numeric_values.append(value)
print(
    "Generated H->gg: 2 orientations, exact color, 10 IBP targets; all tensor UV poles cancel.",
    flush=True,
)
print(
    "IBP:",
    solution.stats,
    "finite coefficient:",
    finite_factor,
    "heavy series:",
    heavy_series,
    flush=True,
)
print(
    "Native sign/tree phase, Ward identities, spin/color norm, width and 19 spacelike/threshold/timelike OneLOop checks passed.",
    flush=True,
)

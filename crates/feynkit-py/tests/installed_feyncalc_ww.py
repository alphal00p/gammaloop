"""Generated massive e- e+ -> W- W+ with Higgs and gauge interference."""

import numpy as np
from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
vertices = [
    v
    for v in model.vertex_rules
    if set(v.particles) <= {"e-", "e+", "a", "Z", "W-", "W+", "ve", "ve~", "H"}
]
generated = model.generate_diagrams(
    [11, -11],
    [-24, 24],
    max_vertices=2,
    maximum_bridges=None,
    vertex_allow=vertices,
    numerator_grouping=None,
    progress=None,
)
assert len(generated.diagrams) == 4
assert {d.internal_edges[0].particle_pdg for d in generated.diagrams} == {
    22,
    23,
    25,
    12,
}
P = S("gammalooprs::P")
charge, sw, cw, mw, mz, mh, me, ye, vev, yme = S(
    "UFO::ee",
    "UFO::sw",
    "UFO::cw",
    "UFO::MW",
    "UFO::MZ",
    "UFO::MH",
    "UFO::Me",
    "UFO::ye",
    "UFO::vev",
    "UFO::yme",
)
s, t, u, H = S("ww::s", "ww::t", "ww::u", "ww::Higgs")
zero, one, pi = E("0"), E("1"), Symbol.PI
kin = hep.Kinematics.mandelstam(
    [P(i) for i in range(4)], [me**2, me**2, mw**2, mw**2], [s, t, u]
)
ports = S("ww::port0", "ww::port1", "ww::port2", "ww::port3")
a, b, c, inverse, wave, rep, index = S(
    "a_", "b_", "c_", "inverse_", "wave_", "rep_", "index_"
)
conjugate, wrapped = S("spenso::conj", "ww::adjoint")
diagram_terms = {}
for diagram in generated.diagrams:
    numerator = model.expand_couplings(
        diagram.numerator_expression(in_lmb=True).to_expression()
    )
    numerator = (
        numerator.replace(ye, model.parameter("ye").expression)
        .replace(yme, me)
        .replace(vev, model.parameter("vev").expression)
        .replace(E("1/2").sqrt(), E("2").sqrt() / 2)
    )
    for edge in diagram.external_edges:
        matches = list(
            diagram.projector_expression().match(
                wave(edge.id, rep(4, index)), max_level=0
            )
        )
        assert len(matches) == 1
        match = dict(matches[0])
        numerator = numerator.replace(
            match[rep](4, match[index]), match[rep](4, ports[edge.external_index])
        )
    denominator = kin.apply(
        diagram.denominator_expression(dimension=4, in_lmb=True)
        .to_expression()
        .replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
    ).expand()
    pdg = diagram.internal_edges[0].particle_pdg
    term = (
        numerator
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
        / denominator
    )
    diagram_terms[pdg] = term

polar_ports = S("ww::physical2", "ww::physical3")
projection = model.particle_by_pdg(-24).spin_sum(
    P(2), ports[2], polar_ports[0]
) * model.particle_by_pdg(24).spin_sum(P(3), ports[3], polar_ports[1])
# Certify that the Z propagator's longitudinal numerator vanishes against
# both physical W sums before removing it. This keeps the later traces small.
longitudinal = (diagram_terms[23] * (s - mz**2)).expand().coefficient(mz**-2) * mz**-2
transverse_test = (longitudinal * projection).replace(
    P(3, a), P(0, a) + P(1, a) - P(2, a)
)
transverse_test = (
    TensorExpression(transverse_test).expand().simplify_metrics().to_dots()
)
transverse_test = (
    kin.apply(transverse_test.to_expression())
    .replace(u, 2 * me**2 + 2 * mw**2 - s - t)
    .expand()
)
assert transverse_test == zero

diagram_terms[23] = (diagram_terms[23] - longitudinal / (s - mz**2)).together().expand()
spins = model.particle_by_pdg(11).spin_sum(
    P(0), ports[0], wrapped(ports[1]), average=True
) * model.particle_by_pdg(-11).spin_sum(P(1), wrapped(ports[0]), ports[1], average=True)
density = model.particle_by_pdg(-24).spin_sum(
    P(2), ports[2], wrapped(ports[2])
) * model.particle_by_pdg(24).spin_sum(P(3), ports[3], wrapped(ports[3]))
operators = {}
adjoints = {}
for pdg, term in diagram_terms.items():
    operator = TensorExpression(term.expand())
    adjoint = (
        operator.dirac_adjoint()
        .expand()
        .simplify_gamma0()
        .to_expression()
        .replace(conjugate(P(a, b)), P(a, b))
    )
    for real in (charge, sw, cw, mw, mz, mh, me, s, t, u, H):
        adjoint = adjoint.replace(conjugate(real), real)
    operators[pdg] = operator
    adjoints[pdg] = TensorExpression(adjoint).wrap_indices(wrapped)
squared = zero
pair_results = {}
for left, left_operator in operators.items():
    for right, right_adjoint in adjoints.items():
        paired = TensorExpression(
            left_operator.to_expression() * right_adjoint * spins,
            cook_indices=CookSettings.indices(),
        )
        paired = paired.expand()
        traced = paired.simplify_gamma().expand().simplify_gamma().simplify_epsilon()
        scalar = (
            TensorExpression(
                traced.to_expression() * density, cook_indices=CookSettings.indices()
            )
            .expand()
            .simplify_metrics()
            .to_dots()
        )
        assert scalar.is_scalar
        result = (
            kin.apply(scalar.to_expression())
            .replace(u, 2 * me**2 + 2 * mw**2 - s - t)
            .replace(mz, mw / cw)
            .replace(cw, (1 - sw**2).sqrt())
            .together()
        )
        # Scalar-product substitutions alone do not impose the dependence of
        # four vectors inside epsilon. Reuse Spenso's compact/indexed conversion.
        conserved = TensorExpression(result).undo_schoonschip().to_expression()
        conserved = conserved.replace(P(3, a), P(0, a) + P(1, a) - P(2, a))
        result = (
            TensorExpression(conserved)
            .expand()
            .simplify_epsilon()
            .simplify_metrics()
            .to_dots()
            .to_expression()
            .together()
        )
        pair_results[left, right] = result
        squared += result * H ** (int(left == 25) + int(right == 25))
    print("Contracted interference row", left, flush=True)
squared = squared.together()

for (left, right), result in pair_results.items():
    assert (result - pair_results[right, left]).together() == zero, (left, right)
print(
    "All four generated diagrams and sixteen interference products contracted",
    flush=True,
)

# Published full massive squared amplitude, including Higgs exchange.
# https://feyncalc.github.io/FeynCalcExamples/EW/Tree/AnelEl-WW
reference = E("""
- ((pi^2*alpha^2*((2*s^2*(s- mh^2)^2*me^8+ 4*s*(s- mh^2)*((- ((s- 4*t*sw^2)*mw^2)-
2*s*t*(sw^2- 1))*mh^2+ s*((- 4*t*sw^2+ s+ 2*t)*mw^2+ s*t*(2*sw^2- 1)))*me^6+
2*(((96*t^2*sw^4- 16*s*t*sw^2- 3*s^2)*mw^4- 2*s*t*(16*t*sw^4+ 4*(2*s- 3*t)*sw^2+
3*s)*mw^2+ s^2*t*(8*t*sw^4- 12*t*sw^2+ s+ 6*t))*mh^4- 2*s*((96*t^2*sw^4- 16*t*(s+
3*t)*sw^2+ s*(2*t- 3*s))*mw^4- s*t*(32*t*sw^4+ 8*(2*s- 3*t)*sw^2+ 9*s+ 4*t)*mw^2+
s^2*t*(8*t*sw^4- 8*t*sw^2+ s+ 4*t))*mh^2+ s^2*((96*t^2*sw^4- 16*t*(s+ 6*t)*sw^2- 3*s^2+
24*t^2+ 4*s*t)*mw^4- 4*s*t*(8*t*sw^4+ (4*s- 6*t)*sw^2+ 3*s+ 4*t)*mw^2+ s^2*t*(8*t*sw^4-
4*t*sw^2+ s+ 4*t)))*me^4- (4*(2*(- 48*t^2*sw^4- 2*s*t*sw^2+ s^2)*mw^6+ 2*t*sw^2*(s*(5*s-
4*t)- 24*(s- 2*t)*t*sw^2)*mw^4+ s*t*(8*(s- 4*t)*t*sw^4+ 12*t^2*sw^2- s*(2*s+ 3*t))*mw^2+
s^2*t^2*(4*(s+ 2*t)*sw^4- 2*(s+ 3*t)*sw^2+ s+ 2*t))*mh^4- 4*s*(4*(- 48*t^2*sw^4- 2*(s-
6*t)*t*sw^2+ s*(s+ t))*mw^6+ 2*t*(- 48*(s- 2*t)*t*sw^4+ 2*(5*s^2- 10*t*s- 12*t^2)*sw^2-
s*(s+ t))*mw^4- s*t*(- 16*(s- 4*t)*t*sw^4+ 4*(s- 6*t)*t*sw^2+ 4*s^2+ 2*t^2+ 5*s*t)*mw^2+
s^2*t^2*(8*(s+ 2*t)*sw^4- 2*(s+ 4*t)*sw^2+ s+ 3*t))*mh^2+ s^2*(8*(- 48*t^2*sw^4- 2*(s-
12*t)*t*sw^2+ s*(s+ 2*t))*mw^6+ 4*t*(- 48*(s- 2*t)*t*sw^4+ 2*(5*s^2- 16*t*s-
24*t^2)*sw^2+ s*(t- 2*s))*mw^4- 4*s*t*(- 8*(s- 4*t)*t*sw^4+ 4*(s- 3*t)*t*sw^2+ 2*s^2+
2*t^2+ 3*s*t)*mw^2+ s^2*t^2*(16*(s+ 2*t)*sw^4- 8*t*sw^2+ s+ 4*t)))*me^2+ 2*(s-
mh^2)^2*(4*(24*t^2*sw^4+ 4*s*t*sw^2+ s^2)*mw^8- 8*t*(4*t*(s+ 6*t)*sw^4+ s*(3*t-
4*s)*sw^2+ s^2)*mw^6+ t*(8*t*(17*s^2+ 20*t*s+ 12*t^2)*sw^4- 20*s^2*t*sw^2+ s^2*(4*s+
5*t))*mw^4- 2*s*t^2*(8*(2*s^2+ 3*t*s+ 2*t^2)*sw^4- 4*(2*s^2+ 2*t*s+ t^2)*sw^2+ s*(2*s+
t))*mw^2+ s^2*t^3*(s+ t)*(8*sw^4- 4*sw^2+ 1)))*mw^4- 2*s*(1- sw^2)*(2*s^2*(s-
mh^2)^2*me^8+ 2*s*(s- mh^2)*(((4*t*sw^2- 2*s+ 2*t)*mw^2+ s*t*(3- 2*sw^2))*mh^2+ s*(2*(-
2*t*sw^2+ s+ t)*mw^2+ s*t*(2*sw^2- 1)))*me^6+ 2*(((s*(2*t- 3*s)- 8*(s-
3*t)*t*sw^2)*mw^4+ s*t*((4*t- 8*s)*sw^2- 5*s+ 6*t)*mw^2+ s^2*t*(- 4*t*sw^2+ s+
3*t))*mh^4+ s*(2*(3*s^2+ 8*t*sw^2*s- 4*t*s+ 6*t^2)*mw^4+ 4*s*t*((4*s- 2*t)*sw^2+ 4*s-
t)*mw^2+ s^2*t*(4*t*sw^2- 2*s- 3*t))*mh^2+ s^2*((- 3*s^2+ 6*t*s+ 12*t^2- 8*t*(s+
3*t)*sw^2)*mw^4- s*t*((8*s- 4*t)*sw^2+ 11*s+ 10*t)*mw^2+ s^2*t*(s+ 2*t)))*me^4+ (-
2*(4*(s*(s+ t)- t*(s+ 12*t)*sw^2)*mw^6+ 2*t*((5*s^2- 16*t*s+ 24*t^2)*sw^2+ s*(5*s+
t))*mw^4- s*t*(4*s^2+ t*s- 6*t^2+ 4*t*(t- s)*sw^2)*mw^2+ s^2*t^2*(- 2*t*sw^2+ s+
t))*mh^4+ s*(8*(2*s^2+ 4*t*s+ 3*t^2- 2*t*(s+ 6*t)*sw^2)*mw^6+ 4*t*(8*s^2- 3*t*s- 6*t^2+
2*(5*s^2- 22*t*s+ 12*t^2)*sw^2)*mw^4- 2*s*t*(8*s^2+ t*s- 8*t^2- 4*(s- 2*t)*t*sw^2)*mw^2+
s^2*t^2*(4*s*sw^2+ s+ 2*t))*mh^2- 4*s^2*(2*(s^2- t*sw^2*s+ 3*t*s+ 3*t^2)*mw^6- t*(-
3*s^2+ (28*t- 5*s)*sw^2*s+ t*s+ 6*t^2)*mw^4- s*t*(2*s^2+ t*s- t^2+ 2*t^2*sw^2)*mw^2+
s^2*t^2*(s+ t)*sw^2))*me^2+ 4*(s- mh^2)^2*mw^2*(2*(2*t*(s+ 3*t)*sw^2+ s*(s+ t))*mw^6-
t*((- 8*s^2+ 10*t*s+ 24*t^2)*sw^2+ 3*s*t)*mw^4+ 2*t*(s^3+ 2*t*(3*s^2+ 5*t*s+
3*t^2)*sw^2)*mw^2- s*t^3*(s+ t)*(2*sw^2- 1)))*mw^2+ (2*s^2*(s- mh^2)^2*me^8+ 4*s*(s-
mh^2)*((s*t- (s- 2*t)*mw^2)*mh^2+ s^2*mw^2)*me^6+ 2*(((- 3*s^2+ 4*t*s+ 12*t^2)*mw^4-
4*s*(s- 2*t)*t*mw^2+ s^2*t*(s+ t))*mh^4- 2*s^2*(- 3*(s- 2*t)*mw^4+ t*(4*t- 7*s)*mw^2+
s^2*t)*mh^2+ s^2*((- 3*s^2+ 8*t*s+ 12*t^2)*mw^4- 2*s*t*(5*s+ 4*t)*mw^2+ s^2*t*(s+
t)))*me^4+ (- ((8*s*(s+ 2*t)*mw^6+ 4*t*(10*s^2+ 13*t*s+ 12*t^2)*mw^4- 4*s*t*(2*s^2+ t*s-
2*t^2)*mw^2+ s^3*t^2)*mh^4)- 8*s*mw^2*(- 2*(s^2+ 3*t*s+ 3*t^2)*mw^4- 3*t*(3*s^2+ 3*t*s+
2*t^2)*mw^2+ s*t*(2*s^2+ t*s- t^2))*mh^2+ 8*s^2*mw^2*(- ((s^2+ 4*t*s+ 6*t^2)*mw^4)-
4*s*t*(s+ t)*mw^2+ s^2*t*(s+ t)))*me^2+ 8*(s- mh^2)^2*mw^4*((s^2+ 2*t*s+ 3*t^2)*mw^4+
2*t*(s^2- 2*t*s- 3*t^2)*mw^2+ t*(s^3+ 3*t*s^2+ 5*t^2*s+ 3*t^3)))*(s-
s*sw^2)^2))/(2*s^2*t^2*(s- mh^2)^2*mw^4*(mw^2- s*cw^2)^2*sw^4))
""")
for _name, _symbol in [
    ("s", s),
    ("t", t),
    ("me", me),
    ("mw", mw),
    ("mh", mh),
    ("sw", sw),
    ("cw", cw),
]:
    reference = reference.replace(S(_name), _symbol)
reference = (
    reference.replace(S("pi"), one)
    .replace(S("alpha"), charge**2 / 4)
    .replace(cw, (1 - sw**2).sqrt())
)
full = squared.replace(H, one)
assert (full - reference).together() == zero
print("PASS massive FeynCalc amplitude", flush=True)
z, angle = S("ww::inverse_s", "ww::angle")
high_energy = (
    squared.replace(t, -s * (1 - angle) / 2)
    .replace(s, 1 / z)
    .series(z, 0, 0)
    .to_expression()
    .expand()
)
assert high_energy.replace(H, one).coefficient(z**-2).together() == zero
assert high_energy.replace(H, one).coefficient(z**-1).together() == zero
print("HIGH ENERGY GROWTH", high_energy.coefficient(z**-1).factor(), flush=True)
assert (
    high_energy.coefficient(z**-1)
    - charge**4 * me**2 * (H - 1) ** 2 / (32 * mw**4 * sw**4)
).together() == zero
massless = full.replace(me, zero).together()
# Integrate the rational t dependence with opaque parameter coefficients,
# then restore those coefficients and differentiate the result exactly.
coefficients = (full * t**2).together().expand().coefficient_list(t)
C = S("ww::C")
compressed = sum(
    (C(i) * power / t**2 for i, (power, _) in enumerate(coefficients)), zero
)
primitive = compressed.integrate(t)
for i, (_, coefficient) in enumerate(coefficients):
    primitive = primitive.replace(C(i), coefficient.together())
assert (primitive.derivative(t) - full).together() == zero
primitive = primitive.replace(t.log(), (-t).log())
print("PASS full massive antiderivative", flush=True)
P = S("gammalooprs::P")
kin = hep.Kinematics.mandelstam(
    [P(i) for i in range(4)], [me**2, me**2, mw**2, mw**2], [s, t, u]
)
flux = kin.flux(P(0), P(1))
phase_space = kin.two_body_phase_space(P(2), P(3))
span = ((s - 4 * me**2) * (s - 4 * mw**2)).sqrt()
# Verify the positive physical branch of (dPhi/dcos)/flux/(dt/dcos).
native_dt_factor = 2 * pi * phase_space / flux / (span / 2)
dt_factor = one / (16 * pi * s * (s - 4 * me**2))
assert (native_dt_factor**2 - dt_factor**2).together() == zero
print("PASS massive flux and phase-space normalization", flush=True)
t_center = me**2 + mw**2 - s / 2
t_upper, t_lower = t_center + span / 2, t_center - span / 2
# Keep the expression in endpoint form; no expansion of the logarithm is needed.
total = dt_factor * (primitive.replace(t, t_upper) - primitive.replace(t, t_lower))
beta = (1 - 4 * mw**2 / s).sqrt()
den = mw**2 - s * (1 - sw**2)
reference_massless = (charge**4 / (16 * pi)) * (
    beta
    * (
        16 * (3 * mw**2 + 8 * s) * (2 * mw**4 + s**2) * sw**2
        - 3 * s * (-20 * s * mw**2 + 32 * mw**4 + 21 * s**2)
        - 4 * (8 * s**2 * mw**2 + 160 * s * mw**4 + 96 * mw**6 + 15 * s**3) * sw**4
    )
    / (96 * s**2 * sw**4 * den**2)
    + ((s - 2 * mw**2 - s * beta) / (s - 2 * mw**2 + s * beta)).log()
    * (
        24 * s * (s * mw**2 + 4 * mw**4 + s**2)
        - 24 * (2 * s**2 * mw**2 + 10 * s * mw**4 + 4 * mw**6 + s**3) * sw**2
    )
    / (96 * s**3 * sw**4 * den)
)
# Independent quadrature uses the published amplitude, not the antiderivative.
_nodes, _weights = np.polynomial.legendre.leggauss(256)
numeric_checks = []
for _electron_mass in (0.0, 0.1, 0.4):
    for _sv in (5.0, 10.0, 25.0):
        _values = {
            charge: 1.0,
            sw: np.sqrt(0.23),
            mw: 1.0,
            mh: 1.5,
            me: _electron_mass,
            s: _sv,
            H: 1.0,
        }
        _native_factor = complex(native_dt_factor.evaluate(_values))
        _physical_factor = complex(dt_factor.evaluate(_values))
        assert (
            _native_factor.real > 0 and abs(_native_factor - _physical_factor) < 1e-12
        )
        _span = np.sqrt((_sv - 4 * _electron_mass**2) * (_sv - 4))
        _center = _electron_mass**2 + 1 - _sv / 2
        _norm = 1 / (16 * np.pi * _sv * (_sv - 4 * _electron_mass**2))
        _quadrature = (
            sum(
                float(w)
                * complex(
                    reference.evaluate({**_values, t: _center + _span * float(x) / 2})
                )
                for x, w in zip(_nodes, _weights)
            )
            * _span
            / 2
            * _norm
        )
        _actual = complex(total.evaluate(_values))
        assert abs(_actual - _quadrature) < 1e-10, (_values, _actual, _quadrature)
        if _electron_mass == 0:
            _published = complex(reference_massless.evaluate(_values))
            assert abs(_actual - _published) < 1e-12, (_actual, _published)
        assert _actual.real > 0 and abs(_actual.imag) < 1e-12
        numeric_checks.append(
            (_electron_mass, _sv, _actual.real, abs(_actual - _quadrature))
        )
print(
    "PASS all nine massive rates and massless total reference",
    numeric_checks,
    flush=True,
)

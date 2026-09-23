"""Generated muon decay, recursive phase space and exact mass correction.

Reference: FeynCalc EW/Tree/Mu-ElAnelNmu. The full unitary-gauge squared
amplitude is retained before the low-energy limit and analytic integrations.
"""

import numpy as np
from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
muon, electron = model.particle_by_pdg(13), model.particle_by_pdg(11)
vertices = [
    v
    for v in model.vertex_rules
    if any(p in ("W-", "W+") for p in v.particles)
    and any(p in ("e-", "e+", "mu-", "mu+") for p in v.particles)
]
assert len(vertices) == 4
generated = model.generate_diagrams(
    [13],
    [11, -12, 14],
    max_vertices=2,
    maximum_bridges=None,
    vertex_allow=vertices,
    numerator_grouping=None,
    progress=None,
)
assert len(generated.diagrams) == 1
diagram = generated.diagrams[0]
assert len(diagram.internal_edges) == 1
assert abs(diagram.internal_edges[0].particle_pdg) == 24
P = S("gammalooprs::P")
charge, sw, mw, mm, me = S("UFO::ee", "UFO::sw", "UFO::MW", "UFO::MM", "UFO::Me")
GF, s, t, z = S(
    "muon_decay::GF", "muon_decay::s", "muon_decay::t", "muon_decay::inverse_W2"
)
M = S("muon_decay::M", is_positive=True)
mass = S("muon_decay::m", is_positive=True)
conjugate, wrapped = S("spenso::conj", "muon_decay::adjoint")
ports = S("muon_decay::mu", "muon_decay::e", "muon_decay::antinue", "muon_decay::numu")
wave, rep, index = S("muon_decay::wave_", "muon_decay::rep_", "muon_decay::index_")
zero, one, pi = E("0"), E("1"), Symbol.PI
numerator = model.expand_couplings(
    diagram.numerator_expression(in_lmb=True).to_expression()
)
for edge in diagram.external_edges:
    matches = list(
        diagram.projector_expression().match(wave(edge.id, rep(4, index)), max_level=0)
    )
    assert len(matches) == 1
    numerator = numerator.replace(dict(matches[0])[index], ports[edge.external_index])
operator = TensorExpression(
    (
        numerator
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    ).expand()
)
adjoint = operator.dirac_adjoint().expand().simplify_gamma0().to_expression()
momentum_index = S("muon_decay::momentum_index_")
adjoint = adjoint.replace(conjugate(P(momentum_index, index)), P(momentum_index, index))
for real in (charge, sw, mw, mm, me):
    adjoint = adjoint.replace(conjugate(real), real)
adjoint = TensorExpression(adjoint).wrap_indices(wrapped)
projector = (
    muon.spin_sum(P(0), ports[0], wrapped(ports[3]), average=True)
    * electron.spin_sum(P(1), wrapped(ports[2]), ports[1])
    * model.particle_by_pdg(-12).spin_sum(P(2), ports[2], wrapped(ports[1]))
    * model.particle_by_pdg(14).spin_sum(P(3), wrapped(ports[0]), ports[3])
)
contracted = (
    TensorExpression(
        operator.to_expression() * adjoint * projector,
        cook_indices=CookSettings.indices(),
    )
    .expand()
    .simplify_gamma()
    .expand()
    .simplify_gamma()
    .simplify_epsilon()
    .expand()
    .simplify_metrics()
    .to_dots()
)
assert contracted.is_scalar
# s=(electron+antinue)^2; t=(electron+numu)^2. All four external momenta
# are future directed, with P0=P1+P2+P3 and both neutrinos massless.
kin = hep.Kinematics()
for i, mass_squared in enumerate((mm**2, me**2, zero, zero)):
    kin = kin.with_scalar_product(P(i), P(i), mass_squared)
for i, j, value in [
    (1, 2, (s - me**2) / 2),
    (1, 3, (t - me**2) / 2),
    (2, 3, (mm**2 + me**2 - s - t) / 2),
    (0, 1, (s + t) / 2),
    (0, 2, (mm**2 - t) / 2),
    (0, 3, (mm**2 - s) / 2),
]:
    kin = kin.with_scalar_product(P(i), P(j), value)
den = S("gammalooprs::denom")
edge_, momentum_, mass_, quad_ = S(
    "muon_decay::edge_",
    "muon_decay::momentum_",
    "muon_decay::mass_",
    "muon_decay::quad_",
)
denominator = kin.apply(
    diagram.denominator_expression(dimension=4, in_lmb=True)
    .to_expression()
    .replace(den(edge_, momentum_, mass_, quad_), quad_)
)
assert (denominator - s + mw**2).expand() == zero
squared = (kin.apply(contracted.to_expression()) / denominator**2).together()
squared = (
    squared.replace(charge**4, 32 * GF**2 * mw**4 * sw**4)
    .replace_multiple([Replacement(mm, M), Replacement(me, mass)])
    .together()
)

# Independent full unitary-gauge reference in scalar products, before taking
# the Fermi limit. This tests the longitudinal W numerator as well as its pole.
pk, pq1, pq2 = (s + t) / 2, (M**2 - t) / 2, (M**2 - s) / 2
kq1, kq2, q1q2 = (s - mass**2) / 2, (t - mass**2) / 2, (M**2 + mass**2 - s - t) / 2
reference_full = (
    16
    * GF**2
    / (s - mw**2) ** 2
    * (
        -2 * mass**2 * pq2 * kq1**2
        - 2 * mass**2 * mw**2 * kq2 * pq1
        + 2 * mass**2 * mw**2 * kq1 * pq2
        - 2 * mass**2 * mw**2 * pk * q1q2
        - mass**4 * kq1 * pq2
        + 2 * mass**2 * pk * kq1 * kq2
        + 2 * mass**2 * kq1 * kq2 * pq1
        + 2 * mass**2 * pk * kq1 * q1q2
        + 2 * mass**2 * kq1 * pq1 * q1q2
        - 4 * mass**2 * mw**2 * pq1 * q1q2
        + 4 * mw**4 * kq2 * pq1
    )
)
assert (squared - reference_full).together() == zero
fermi = squared.replace(mw, 1 / z.sqrt()).series(z, 0, 0).to_expression().expand()
assert (fermi - 16 * GF**2 * (t - mass**2) * (M**2 - t)).together() == zero
finite_W = (
    squared.replace(mass, zero)
    .replace(mw, 1 / z.sqrt())
    .series(z, 0, 1)
    .to_expression()
    .expand()
)
assert (finite_W - fermi.replace(mass, zero) * (1 + 2 * s * z)).together() == zero

# Factor P -> Q+antinue, Q -> electron+numu with Q^2=t. The shared two-body
# measures retain their angles; integrate the global orientation (4pi) and
# the inner azimuth (2pi), and include the recursive dt/(2pi).
Q = S("muon_decay::Q")
outer = (
    hep.Kinematics()
    .with_scalar_product(Q, Q, t)
    .with_scalar_product(P(2), P(2), zero)
    .with_scalar_product(Q, P(2), (M**2 - t) / 2)
)
inner = (
    hep.Kinematics()
    .with_scalar_product(P(1), P(1), mass**2)
    .with_scalar_product(P(3), P(3), zero)
    .with_scalar_product(P(1), P(3), (t - mass**2) / 2)
)
outer_measure = outer.two_body_phase_space(Q, P(2))
inner_measure = inner.two_body_phase_space(P(1), P(3))
# Select the positive roots on mass^2 < t < M^2.
outer_measure = outer_measure.replace(
    (4 * ((M**2 - t) / 2).expand() ** 2).sqrt(), M**2 - t
).together()
inner_measure = inner_measure.replace(
    (4 * ((t - mass**2) / 2).expand() ** 2).sqrt(), t - mass**2
).together()
assert (outer_measure - (M**2 - t) / (32 * pi**2 * M**2)).together() == zero
assert (inner_measure - (t - mass**2) / (32 * pi**2 * t)).together() == zero
cosine = S("muon_decay::cosine")
electron_energy = (t + mass**2) / (2 * t.sqrt())
electron_momentum = (t - mass**2) / (2 * t.sqrt())
spectator_energy = (M**2 - t) / (2 * t.sqrt())
s_of_cosine = mass**2 + 2 * spectator_energy * (
    electron_energy - electron_momentum * cosine
)
s_min = s_of_cosine.replace(cosine, one).together()
s_max = s_of_cosine.replace(cosine, -one).together()
assert (s_min - mass**2 * M**2 / t).together() == zero
assert (s_max - M**2 - mass**2 + t).together() == zero
jacobian = -s_of_cosine.derivative(cosine)
recursive_measure = (
    outer_measure * inner_measure * (4 * pi) * (2 * pi) / (2 * pi * jacobian)
).together()
dalitz_measure = kin.three_body_phase_space(P(1), P(2), P(3)).replace(mm, M).together()
assert (dalitz_measure - recursive_measure).together() == zero
assert (dalitz_measure - 1 / (128 * pi**3 * M**2)).together() == zero
flux = hep.Kinematics().with_scalar_product(P(0), P(0), M**2).flux(P(0))
assert flux == 2 * M
density = (fermi * dalitz_measure / flux).together()

# Exact finite-electron-mass correction, using dimensionless t/M^2 and m^2/M^2.
tau = S("muon_decay::tau", is_positive=True)
r = S("muon_decay::r", is_positive=True)
normalization = GF**2 * M**5 / (192 * pi**3)
integrand = (
    (density * (s_max - s_min) * M**2 / normalization)
    .replace(t, M**2 * tau)
    .replace(mass, M * r.sqrt())
    .together()
)
assert (integrand - 12 * (1 - tau) ** 2 * (tau - r) ** 2 / tau).together() == zero
primitive = integrand.integrate(tau)
mass_correction = (primitive.replace(tau, one) - primitive.replace(tau, r)).expand()
reference_correction = 1 - 8 * r + 8 * r**3 - r**4 - 12 * r**2 * r.log()
assert (mass_correction - reference_correction).together() == zero
assert mass_correction.replace(r, one) == zero
total_width = normalization * mass_correction

# Derive the Michel spectrum with xE=2*E_e/M=(s+t)/M^2 in the massless limit.
xE = S("muon_decay::xE", is_positive=True)
massless_density = density.replace(mass, zero)
dimensionless_density = (
    (massless_density * M**4 / normalization).replace(t, M**2 * tau).together()
)
michel_primitive = dimensionless_density.integrate(tau)
michel_spectrum = (
    michel_primitive.replace(tau, xE) - michel_primitive.replace(tau, zero)
).expand()
assert (michel_spectrum - 2 * xE**2 * (3 - 2 * xE)).together() == zero
michel_integral = michel_spectrum.integrate(xE)
assert (
    michel_integral.replace(xE, one) - michel_integral.replace(xE, zero) - 1
).together() == zero

# Integrate the first finite-W correction from the generated full amplitude.
correction_density = (finite_W.coefficient(z) * dalitz_measure / flux).replace(
    mass, zero
)
correction_primitive = correction_density.integrate(s)
correction_t = (
    correction_primitive.replace(s, M**2 - t) - correction_primitive.replace(s, zero)
).together()
correction_primitive = correction_t.integrate(t)
W_coefficient = (
    correction_primitive.replace(t, M**2) - correction_primitive.replace(t, zero)
).together()
assert (W_coefficient / normalization - 3 * M**2 / 5).together() == zero

# Independent quadrature on the physical Dalitz range, including tau-like
# daughter masses; transformed intervals avoid evaluating log(0).
nodes, weights = np.polynomial.legendre.leggauss(160)
numeric_checks = []
for rv in [0.0, (0.511 / 105.658) ** 2, 0.01, 0.1, 0.3, 0.7]:
    value = (
        1.0 if rv == 0 else complex(mass_correction.replace(r, E(str(rv))).evaluate({}))
    )
    assert abs(complex(value).imag) < 1e-14
    actual = complex(value).real
    integral_value = sum(
        float(w)
        * (1 - rv)
        / 2
        * 12
        * (1 - (rv + (1 - rv) * (float(n) + 1) / 2)) ** 2
        * ((rv + (1 - rv) * (float(n) + 1) / 2) - rv) ** 2
        / (rv + (1 - rv) * (float(n) + 1) / 2)
        for n, w in zip(nodes, weights, strict=True)
    )
    assert abs(actual - integral_value) < 2e-8, (rv, actual, integral_value)
    assert 0 < actual <= 1
    numeric_checks.append((rv, actual, integral_value))
print(
    "Generated full muon amplitude, recursive phase space, Michel spectrum, exact mass correction and finite-W term passed"
)

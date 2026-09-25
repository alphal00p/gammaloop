"""Generated massive photon self-energy through shared tensor and OneLOop APIs.

Run in the installed symbolica.community.hep host (which includes OneLOop).
The photon UV pole is a component of the FeynCalc QED renormalization example:
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization
The full renormalization example, including its other diagrams, remains pending.
"""

from math import pi, sqrt
from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
generated = (
    fk.Process(model, [22], [22])
    .with_loop_count(1, 1)
    .generate_diagrams(
        max_vertices=2,
        maximum_bridges=None,
        vertex_allow=["V_98"],
        numerator_grouping=None,
        progress=None,
    )
)
assert len(generated.diagrams) == 1
diagram = generated.diagrams[0]
D, s, mass, charge, mu, nu = S("D", "s", "UFO::Me", "UFO::ee", "mu", "nu")
K, P, mink, metric = S("gammalooprs::K", "gammalooprs::P", "spenso::mink", "spenso::g")
index, wave = S("index_", "wave_")
numerator = model.expand_couplings(
    diagram.numerator_expression(in_lmb=True).to_expression()
)
for edge in diagram.external_edges:
    matches = list(
        diagram.projector_expression().match(wave(edge.id, mink(4, index)), max_level=0)
    )
    assert len(matches) == 1
    port = dict(matches[0])[index]
    numerator = numerator.replace(mink(4, port), mink(D, (mu, nu)[edge.external_index]))
# Continue Lorentz contractions to D before taking the Dirac trace; tr(1) = 4.
numerator = numerator.replace(mink(4, index), mink(D, index))
trace = TensorExpression(numerator.expand()).simplify_gamma().expand().to_dots()
kinematics = fk.Kinematics(D, momenta=[K(0), P(0)]).with_scalar_product(P(0), P(0), s)
family = diagram.integral_family(kinematics=kinematics)
reducer = (
    fk.TensorReducer(D)
    .with_integrated_vector(K(0, mink(D)))
    .with_external_vector(P(0, mink(D)))
)
reduced = kinematics.apply(reducer.reduce(trace.to_expression()))
d0, d1, a0, b0 = S("d0", "d1", "a0", "b0")
scalar = (
    family.rewrite_numerator(reduced, [d0, d1])
    * diagram.overall_factor_expression(evaluate=True)
    * diagram.numerator_prefactor_expression()
    / (d0 * d1)
)

# The equal-mass bubble needs only A0 and B0. Prove the shifted tadpole
# moments with the existing family map and vacuum projector, before using them.
shift = family.sector([0, 1]).mapping_to(family.sector([1, 0]), [K(0) + P(0)])
assert shift is not None
vacuum = fk.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
for moment in (shift.apply(family.denominators[0]), family.denominators[1]):
    averaged = kinematics.apply(vacuum.reduce(moment))
    assert (averaged - family.denominators[0] - s).expand() == 0
weights = {
    E("1"): E("0"),  # Scaleless polynomial moment in dimensional regularization.
    1 / d0: a0,
    1 / d1: a0,
    d0 / d1: s * a0,
    d1 / d0: s * a0,
    1 / (d0 * d1): b0,
}
parts = scalar.expand().coefficient_list(d0, d1)
assert {monomial for monomial, _ in parts} <= weights.keys()
integrated = sum(coefficient * weights[monomial] for monomial, coefficient in parts)
transverse = s * metric(mink(D, mu), mink(D, nu)) - P(0, mink(D, mu)) * P(
    0, mink(D, nu)
)
form_factor = (
    integrated.expand().coefficient(metric(mink(D, mu), mink(D, nu))) / s
).together()
expected = (
    2
    * charge**2
    * (2 * (D - 2) * a0 - ((D - 2) * s + 4 * mass**2) * b0)
    / (s * (D - 1))
)
assert (form_factor - expected).together() == 0
assert (integrated - form_factor * transverse).together() == 0

# OneLOop coefficient functions remain ordinary evaluable Symbolica expressions.
# Positive mass and scale; s != 0 because the external Gram matrix is inverted.
eps, a_finite, b_finite, scale2 = S("eps", "a_finite", "b_finite", "scale2")
a_master = S("oneloopmaster::A0")(mass**2, scale2)
b_master = S("oneloopmaster::B0")(s, mass**2, mass**2, scale2)
a_coefficients = oneloop.master_coefficients(a_master)
b_coefficients = oneloop.master_coefficients(b_master)
laurent = (
    form_factor.replace(D, 4 - 2 * eps)
    .replace(a0, a_finite + mass**2 / eps)
    .replace(b0, b_finite + 1 / eps)
    .series(eps, 0, 0)
    .to_expression()
)
assert (laurent.coefficient(eps**-1) + 4 * charge**2 / 3).together() == 0
finite = dict(laurent.expand().coefficient_list(eps))[E("1")]
expected_finite = (
    4 * charge**2 / 9
    + 8 * charge**2 * (a_finite - mass**2) / (3 * s)
    - 4 * charge**2 * (1 + 2 * mass**2 / s) * b_finite / 3
)
assert (finite - expected_finite).together() == 0
finite = finite.replace(a_finite, a_coefficients[0]).replace(
    b_finite, b_coefficients[0]
)

# Independent 60-digit quadrature of 8 integral_0^1 dx x(1-x)
# log((m^2-s*x*(1-x)-i0)/mu^2), split at the roots above threshold.
# Only the reference constants are stored; no numerical integration dependency.
references = [
    (-1, 1, 1, 0.24171817135108418518155622477374160828 + 0j),
    (3, 1, 1, -1.31288983076412170282358776645606558180 + 0j),
    (
        5,
        1,
        1,
        -2.48545886575608134962861352071442542810
        - 2.62259749958853785344504653346316072513j,
    ),
    (-3, 2, 4, -0.57706419792394813695321259409264475428 + 0j),
    (
        12,
        2,
        3,
        -2.30000503960920344702651658972063293907
        - 3.22453220308305395661169468025272130184j,
    ),
]
for invariant, mass2, scale, reference in references:
    point = {s: float(invariant), mass: sqrt(mass2), scale2: float(scale), charge: 1.0}
    # Check both scalar pole coefficients used in the exact Laurent expansion.
    assert abs(a_coefficients[1].evaluate(point) - mass2) < 1e-12
    assert abs(b_coefficients[1].evaluate(point) - 1) < 1e-12
    assert abs(a_coefficients[2].evaluate(point)) < 1e-12
    assert abs(b_coefficients[2].evaluate(point)) < 1e-12
    value = finite.evaluate(point)
    assert abs(value - reference) < 3e-12, (point, value, reference)
    cut = (
        -4 * pi / 3 * sqrt(1 - 4 * mass2 / invariant) * (1 + 2 * mass2 / invariant)
        if invariant > 4 * mass2
        else 0
    )
    assert abs(value.imag - cut) < 3e-12
print(
    "Photon self-energy: transverse, exact UV pole, five finite-part and cut checks passed"
)

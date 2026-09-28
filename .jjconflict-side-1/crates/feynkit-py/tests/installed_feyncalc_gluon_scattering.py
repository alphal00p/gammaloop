"""Generated gg -> gg with symbolic D and SU(N), plus four-dimensional cut rates.

Reference: https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/GlGl-GlGl
Spin/color completeness, diagram generation and tensor algebra use shared owners.
The gallery uses fixed 1/2 incoming spin factors; dimensional averages use 1/(D-2).
"""

import numpy as np
from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

index_scope = S("spenso::index_scope")

model = hep.Model.standard_model()
P = S("gammalooprs::P")
s = S("gg::s", is_positive=True)
t, u, D, N, dA, gs = S("gg::t", "gg::u", "gg::D", "gg::N", "gg::dA", "UFO::G")
coad, conj = S("spenso::coad", "spenso::conj")
index, a, b, c, inv, wave, rep = S("index_", "a_", "b_", "c_", "inv_", "wave_", "rep_")
ports = S("gg::i0", "gg::i1", "gg::i2", "gg::i3")
wrapped = S("gg::adjoint")
bars = S("gg::j0", "gg::j1", "gg::j2", "gg::j3")
gluon = model.particle("g")
vertices = [
    v
    for v in model.vertex_rules
    if v.particles.count(gluon.name) == len(v.particles) and len(v.particles) in (3, 4)
]
assert len(vertices) == 2
generated = model.process(
    ["g", "g"], ["g", "g"], vertex_allow=vertices
).generate_diagrams(
    max_vertices=2, maximum_bridges=None, numerator_grouping=None, progress=None
)
assert len(generated.diagrams) == 4
kin = (
    hep.Kinematics(D)
    .with_scalar_product(P(0), P(0), E("0"))
    .with_scalar_product(P(1), P(1), E("0"))
    .with_scalar_product(P(2), P(2), E("0"))
    .with_scalar_product(P(3), P(3), E("0"))
    .with_scalar_product(P(0), P(1), s / 2)
    .with_scalar_product(P(2), P(3), s / 2)
    .with_scalar_product(P(0), P(2), -t / 2)
    .with_scalar_product(P(1), P(3), -t / 2)
    .with_scalar_product(P(0), P(3), -u / 2)
    .with_scalar_product(P(1), P(2), -u / 2)
)
terms = []
denominators = []
for diagram in generated.diagrams:
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
    numerator = TensorExpression(numerator).with_lorentz_dimension(D).to_expression()
    denominator = kin.apply(
        diagram.denominator_expression(dimension=D, in_lmb=True)
        .to_expression()
        .replace(S("gammalooprs::denom")(a, b, c, inv), inv)
    ).expand()
    denominators.append(denominator)
    terms.append(
        numerator
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
        / denominator
    )
assert set(denominators) == {E("1"), s, t, u}
amplitude = sum(terms, E("0"))
operator = kin.apply(TensorExpression(amplitude).contract().to_expression().to_dots())
amplitude = operator.to_expression()
assert len(operator.structure.slots) == 8
adjoint = operator.spenso_conjugate().to_expression().replace(conj(P(a, b)), P(a, b))
for real in (s, t, u, D, gs):
    adjoint = adjoint.replace(conj(real), real)
adjoint = TensorExpression(adjoint).wrap_indices(wrapped).to_expression()
for i in range(4):
    adjoint = adjoint.replace(index_scope(wrapped, ports[i]), bars[i])
colors = E("1")
for i in range(4):
    colors *= gluon.color_sum(ports[i], bars[i], average=i < 2)
generic = (amplitude * adjoint * colors * gluon.color**2 / dA**2).replace(
    coad(8, index), coad(dA, index)
)
print("color contraction", flush=True)
colored = (
    TensorExpression(generic, cook_indices=CookSettings.indices())
    .simplify_color()
    .to_expression()
    .to_expression()
    .replace(dA, N**2 - 1)
)
colored = TensorExpression(colored).to_cof_dimension_invariants().to_expression()
# Physical axial references pair the two incoming and the two outgoing gluons.
result = colored
for i, j in enumerate((1, 0, 3, 2)):
    projector = gluon.spin_sum(P(i), ports[i], bars[i], reference=P(j), dimension=D) / (
        2 if i < 2 else 1
    )
    projector = kin.apply(projector)
    print("polarization", i, "start", flush=True)
    tensor = TensorExpression(result * projector)
    tensor = tensor.contract().to_expression().to_dots()
    result = tensor.to_expression()
    result = kin.apply(result).replace(u, -s - t)
    print("polarization", i, "done", flush=True)
expected = (
    (D - 2) ** 2
    * N**2
    * gs**4
    * (t**2 + t * u + u**2) ** 3
    / ((N**2 - 1) * s**2 * t**2 * u**2)
)
delta = (result - expected.replace(u, -s - t)).together()
assert delta == 0

fixed_two = result.together()
dimensional_average = (4 * fixed_two / (D - 2) ** 2).together()
assert dimensional_average.derivative(D).together() == 0
standard = E("9/2") * gs**4 * (3 - t * u / s**2 - s * u / t**2 - s * t / u**2)
assert (
    fixed_two.replace(D, 4).replace(N, 3) - standard.replace(u, -s - t)
).together() == 0
assert (fixed_two - fixed_two.replace(t, -s - t)).together() == 0
print(
    "Generated four-gluon amplitude: symbolic D, SU(N), both averages and Bose exchange passed",
    flush=True,
)

z, cutoff, alpha = S("gg::z", "gg::cutoff", "gg::alpha_s")
kin4 = hep.Kinematics.mandelstam([P(i) for i in range(4)], [E("0")] * 4, [s, t, u])
phase = kin4.two_body_phase_space(P(2), P(3)) / kin4.flux(P(0), P(1))
assert (phase - 1 / (64 * Symbol.PI**2 * s)).together() == 0
# Azimuth integration supplies 2*pi and identical gluon events supply 1/2!.
# Display s*sigma/alpha_s**2, a dimensionless angular-cut rate.
density = (
    (fixed_two.replace(D, 4) * phase * Symbol.PI * s / alpha**2)
    .replace(gs**4, (4 * Symbol.PI * alpha) ** 2)
    .replace(t, -s * (1 - z) / 2)
    .together()
)
assert density.derivative(s).together() == 0
primitive = density.integrate(z)
assert (primitive.derivative(z) - density).together() == 0
primitive = primitive.replace((z - 1).log(), (1 - z).log())
cut_rate = primitive.replace(z, cutoff) - primitive.replace(z, -cutoff)
assert cut_rate.replace(cutoff, 0).together() == 0
assert (
    cut_rate.derivative(cutoff)
    - density.replace(z, cutoff)
    - density.replace(z, -cutoff)
).together() == 0
nodes, weights = np.polynomial.legendre.leggauss(128)
rate_checks = []
for colors in (2, 3, 5):
    for cut in (0.2, 0.5, 0.8):
        exact = cut_rate.evaluate({N: colors, cutoff: cut})
        numeric = cut * sum(
            weight * density.evaluate({N: colors, z: cut * node})
            for node, weight in zip(nodes, weights, strict=True)
        )
        assert exact.real > 0 and abs(exact.imag) < 1e-10
        assert abs(exact - numeric) < 2e-11 * max(1, abs(exact))
        rate_checks.append((colors, cut, exact.real, abs(exact - numeric)))
print(
    "Native flux, identical-event factor and nine angular-cut quadratures passed",
    flush=True,
)

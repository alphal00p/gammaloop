"""Generated diphoton annihilation, gauge checks and identical-photon phase space.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-GaGa
The ordinary amplitudes retain both labeled diagrams and their interference.
"""

from pathlib import Path

from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, Representation, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
s = S("diphoton::s", is_positive=True)
t, u, mass, charge = S("diphoton::t", "diphoton::u", "UFO::Me", "UFO::ee")
kin = fk.Kinematics.mandelstam(
    [P(0), P(1), P(2), P(3)], [mass**2, mass**2, E("0"), E("0")], [s, t, u]
)
generated = (
    fk.Process(model, [11, -11], [22, 22])
    .with_loop_count(0, 0)
    .generate_diagrams(
        max_vertices=2,
        maximum_bridges=None,
        vertex_allow=["V_98"],
        numerator_grouping=None,
        progress=None,
    )
)
assert len(generated.diagrams) == 2
ports = S(
    "diphoton::external_0",
    "diphoton::external_1",
    "diphoton::external_2",
    "diphoton::external_3",
)
a, b, c, inverse, wave, rep, index = S(
    "a_", "b_", "c_", "inverse_", "wave_", "rep_", "index_"
)
conjugate, adjoint_index = S("spenso::conj", "diphoton::adjoint_index")
amplitude = E("0")
denominators = []
for diagram in generated.diagrams:
    numerator = model.expand_couplings(
        diagram.numerator_expression(in_lmb=True).to_expression()
    )
    # Align external tensor ports with the generated wavefunctions, independent
    # of the two diagrams' internal half-edge numbering.
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
    denominators.append(denominator)
    amplitude += (
        numerator
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
        / denominator
    )
assert set(denominators) == {t - mass**2, u - mass**2}
operator = TensorExpression(amplitude.expand())
assert len(operator.structure.slots) == 4
adjoint = operator.dirac_adjoint().expand().simplify_gamma0().to_expression()
adjoint = adjoint.replace(conjugate(P(a, b)), P(a, b))
# Physical momenta and tree-level parameters are real; keep that assumption explicit.
for real in (mass, charge, s, t, u):
    adjoint = adjoint.replace(conjugate(real), real)
adjoint = TensorExpression(adjoint).wrap_indices(adjoint_index)
# Dirac adjunction exchanges the open fermion chain's endpoints.
spins = model.particle_by_pdg(11).spin_sum(
    P(0), ports[0], adjoint_index(ports[1]), average=True
) * model.particle_by_pdg(-11).spin_sum(
    P(1), adjoint_index(ports[0]), ports[1], average=True
)
spin_summed = (
    TensorExpression(
        operator.to_expression() * adjoint * spins, cook_indices=CookSettings.indices()
    )
    .expand()
    .simplify_gamma()
    .expand()
    .simplify_gamma()
)
x, y = t - mass**2, u - mass**2
expected = (
    2
    * charge**4
    * (x / y + y / x + 4 * mass**2 * s / (x * y) - 4 * mass**4 * s**2 / (x * x * y * y))
)
results = []
# Compare covariant, null-reference and timelike-reference photon sums, then
# replace either photon polarization by its momentum (each Ward identity).
for references, ward in (
    ((None, None), None),
    ((P(3), P(2)), None),
    ((P(0), P(0)), None),
    ((None, None), 2),
    ((None, None), 3),
):
    density = E("1")
    for position, reference in zip((2, 3), references):
        if position == ward:
            mink = Representation.mink(4)
            longitudinal = P(position, mink(ports[position]).to_expression())
            density *= longitudinal * longitudinal.replace(
                ports[position], adjoint_index(ports[position])
            )
        else:
            density *= model.particle_by_pdg(22).spin_sum(
                P(position),
                ports[position],
                adjoint_index(ports[position]),
                reference=reference,
            )
    scalar = (
        TensorExpression(
            spin_summed.to_expression() * density, cook_indices=CookSettings.indices()
        )
        .expand()
        .simplify_metrics()
        .to_dots()
    )
    assert scalar.is_scalar
    squared = (
        kin.apply(scalar.to_expression()).replace(s, 2 * mass**2 - t - u).together()
    )
    target = E("0") if ward else expected.replace(s, 2 * mass**2 - t - u)
    assert (squared - target).together() == E("0"), (references, ward, squared)
    if ward is None:
        results.append(squared)
    print("PASS", references, ward, flush=True)
assert all(result == results[0] for result in results)
assert (
    results[0] - results[0].replace_multiple([Replacement(t, u), Replacement(u, t)])
).together() == E("0")
massless = results[0].replace(mass, E("0"))
assert (massless - 2 * charge**4 * (t / u + u / t)).together() == E("0")
cos_theta, alpha, cutoff = S(
    "diphoton::cos_theta", "diphoton::alpha", "diphoton::cutoff"
)
kin_massless = fk.Kinematics.mandelstam(
    [P(0), P(1), P(2), P(3)], [E("0")] * 4, [s, t, u]
)
angular = (
    massless.replace(t, -s * (1 - cos_theta) / 2)
    .replace(u, -s * (1 + cos_theta) / 2)
    .replace(charge**4, (4 * Symbol.PI * alpha) ** 2)
)
labeled = (
    angular
    * kin_massless.two_body_phase_space(P(2), P(3))
    / kin_massless.flux(P(0), P(1))
).together()
assert (
    labeled - alpha**2 / s * (1 + cos_theta**2) / (1 - cos_theta**2)
).together() == E("0")
# The phase-space API excludes identical-particle factors. Count each event
# once by selecting its forward photon, or equivalently use 1/2! on the full sphere.
primitive = (labeled * s / alpha**2).integrate(cos_theta).together()
assert (primitive.derivative(cos_theta) - labeled * s / alpha**2).together() == E("0")
assert primitive.replace(cos_theta, E("0")) == E("0")
event_cross_section = (
    2 * Symbol.PI * alpha**2 / s * primitive.replace(cos_theta, cutoff)
).together()
assert (
    event_cross_section - 2 * Symbol.PI * alpha**2 / s * (2 * cutoff.atanh() - cutoff)
).together() == E("0")
full_sphere = (
    primitive.replace(cos_theta, cutoff) - primitive.replace(cos_theta, -cutoff)
) / 2
logarithmic = ((1 + cutoff) / (1 - cutoff)).log() - cutoff
for value in ("1/5", "1/2", "4/5"):
    for result in (full_sphere, primitive.replace(cos_theta, cutoff)):
        error = (result - logarithmic).replace(cutoff, E(value)).to_float(40)
        assert abs(complex(error)) < 1e-30
print(
    "Massless angular distribution and identical-photon angular-cut cross section passed"
)

# Retain the electron mass: beta is the incoming CM speed, 0 < beta < 1.
beta = S("diphoton::beta", is_positive=True)
cm_squared = (
    results[0]
    .replace(t, mass**2 - s * (1 - beta * cos_theta) / 2)
    .replace(u, mass**2 - s * (1 + beta * cos_theta) / 2)
    .replace(mass, (s * (1 - beta**2) / 4).sqrt())
).together()
kernel = (cm_squared / (4 * charge**4)).together()
expected_kernel = (1 + beta**2 * cos_theta**2) / (
    1 - beta**2 * cos_theta**2
) + 2 * beta**2 * (1 - beta**2) * (1 - cos_theta**2) / (1 - beta**2 * cos_theta**2) ** 2
assert (kernel - expected_kernel).together() == E("0")
assert kernel.replace(beta, E("0")).together() == E("1")
massive_primitive = kernel.integrate(cos_theta).together()
assert (massive_primitive.derivative(cos_theta) - kernel).together() == E("0")
assert massive_primitive.replace(cos_theta, E("0")).together() == E("0")
expected_primitive = (
    (3 - beta**4) / beta * (beta * cos_theta).atanh()
    - cos_theta
    - (1 - beta**2) ** 2 * cos_theta / (1 - beta**2 * cos_theta**2)
)
assert (massive_primitive - expected_primitive).together() == E("0")
# The massive flux contributes 1/beta; the final photons have unit speed.
# Select its physical positive branch, s > 0 and beta > 0. The shared
# phase-space measure excludes the identical-particle factor.
massive_flux = (
    kin.flux(P(0), P(1)).replace(mass, (s * (1 - beta**2) / 4).sqrt()).expand()
)
assert (massive_flux**2 - 4 * s**2 * beta**2).together() == E("0")
massive_flux = massive_flux.replace((s**2 * beta**2).sqrt(), s * beta)
massive_density_factor = (
    (4 * charge**4 * kin.two_body_phase_space(P(2), P(3)) / massive_flux)
    .replace(charge**4, (4 * Symbol.PI * alpha) ** 2)
    .together()
)
assert (massive_density_factor - alpha**2 / (s * beta)).together() == E("0")
massive_cross_section = (
    2
    * Symbol.PI
    * massive_density_factor
    * massive_primitive.replace(cos_theta, cutoff)
)
assert (
    massive_cross_section.replace(beta, E("1")) - event_cross_section
).together() == E("0")
print(
    "Massive CM distribution, threshold normalization and Symbolica angular integration passed"
)

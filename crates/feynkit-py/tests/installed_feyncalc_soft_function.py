"""Generated light-cone soft radiation with open spinor and photon indices.

Reference:
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/Ga-QQbar-SoftFunction
The reference's explicit 1/(32 pi^2) prefactor is included; this is not an
integrated or dimensionally regulated real-emission phase-space calculation.
"""

from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import (
    GammaChainOrdering,
    GammaSimplifySettings,
    Representation,
    SchoonschipSettings,
    TensorExpression,
    TensorName,
)

model = hep.Model.standard_model()
vertices = [
    vertex
    for vertex in model.vertex_rules
    if sorted(vertex.particles) in [sorted(["b", "b~", "a"]), sorted(["b", "b~", "g"])]
]
assert len(vertices) == 2
P = S("gammalooprs::P")
D, Q, Qbar, kp, km, lam = S(
    "soft::D", "soft::Q", "soft::Qbar", "soft::kp", "soft::km", "soft::lam"
)
# Rank-one declarations let Spenso recognize contracted vectors and slashes.
n, nb, k, ref = [
    TensorName.vector("soft::" + name).to_expression() for name in ("n", "nb", "k", "r")
]
mink, bis = S("spenso::mink", "spenso::bis")
zero, one = E("0"), E("1")
kin = hep.Kinematics(D, momenta=[n, nb, k])
for vector in (n, nb, k):
    kin = kin.with_scalar_product(vector, vector, zero)
kin = (
    kin.with_scalar_product(n, nb, E("2"))
    .with_scalar_product(n, k, kp)
    .with_scalar_product(nb, k, km)
)
ports = S("soft::i0", "soft::i1", "soft::i2", "soft::i3")
a, b, c, inv, idx = S("a_", "b_", "c_", "inv_", "idx_")
gamma, metric = S("spenso::gamma", "spenso::g")
si, sj, sx, sy, sm, sn = S(
    "soft::si", "soft::sj", "soft::sx", "soft::sy", "soft::sm", "soft::sn"
)
ordering = GammaSimplifySettings(chain_ordering=GammaChainOrdering.Canonical)


def simplify(expression):
    return (
        kin.apply(
            TensorExpression(expression)
            .expand()
            .schoonschip(SchoonschipSettings(simplify_chain_like_functions=True))
            .simplify_gamma(ordering)
            .expand()
            .simplify_metrics()
            .to_dots()
        )
        .to_expression()
        .expand()
    )


def projector(left, right, middle):
    # P_minus = nbar_slash n_slash / 4; all spinor dimensions remain four.
    return (
        gamma(bis(4, left), bis(4, middle), nb(mink(D)))
        * gamma(bis(4, middle), bis(4, right), n(mink(D)))
        / 4
    )


def project(expression):
    # The collinear bra and anticollinear ket both select P_minus.
    expression = expression.replace(bis(4, ports[1]), bis(4, sx)).replace(
        bis(4, ports[2]), bis(4, sy)
    )
    return projector(ports[1], sx, sm) * expression * projector(sy, ports[2], sn)


# Check the large-component projector before applying it to generated amplitudes.
pm = projector(si, sj, sm)
identity = metric(bis(4, si), bis(4, sj))
pp = identity - pm
assert simplify(pm.replace(bis(4, sj), bis(4, sx)) * projector(sx, sj, sn) - pm) == zero
assert (
    simplify(
        pm.replace(bis(4, sj), bis(4, sx))
        * pp.replace(bis(4, si), bis(4, sx)).replace(bis(4, sm), bis(4, sn))
    )
    == zero
)
assert simplify(pm.replace(bis(4, sj), bis(4, si))) == 2
assert (
    simplify(
        pm.replace(bis(4, sj), bis(4, sx)) * gamma(bis(4, sx), bis(4, sj), n(mink(D)))
    )
    == zero
)
assert (
    simplify(
        gamma(bis(4, si), bis(4, sx), nb(mink(D))) * pm.replace(bis(4, si), bis(4, sx))
    )
    == zero
)

# The transverse metric projects onto the D-2 directions orthogonal to n,nbar.
mu, nu, rho = S("soft::mu", "soft::nu", "soft::rho")
transverse = (
    metric(mink(D, mu), mink(D, nu))
    - (n(mink(D, mu)) * nb(mink(D, nu)) + nb(mink(D, mu)) * n(mink(D, nu))) / 2
)
assert simplify(transverse * n(mink(D, nu))) == zero
assert simplify(transverse * nb(mink(D, nu))) == zero
assert simplify(transverse * metric(mink(D, mu), mink(D, nu))) == D - 2
assert (
    simplify(transverse * transverse.replace(mu, rho) - transverse.replace(nu, rho))
    == zero
)
print("PASS: symbolic-D transverse metric and large-component projectors")

channels, denominators, generated_channels = [], [], []
for outgoing in (["b", "b~"], ["b", "b~", "g"]):
    result = model.process(["a"], outgoing, vertex_allow=vertices).generate_diagrams(
        max_vertices=len(outgoing) - 1,
        maximum_bridges=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(result.diagrams) == len(outgoing) - 1
    generated_channels.append(result)
    for diagram in result.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        ).replace(S("UFO::MB"), zero)
        for half in diagram.half_edges:
            edge = half.edge.data
            if edge.is_external:
                numerator = numerator.replace(
                    S("gammalooprs::hedge")(half.data, 1), ports[edge.external_index]
                )
        numerator = (
            TensorExpression(numerator).with_lorentz_dimension(D).to_expression()
        )
        denominator = (
            diagram.denominator_expression(dimension=D, in_lmb=True)
            .to_expression()
            .replace(S("gammalooprs::denom")(a, b, c, inv), inv)
            .replace(S("UFO::MB"), zero)
        )
        # Expand the dot products before substituting linear vector combinations.
        denominator = TensorExpression(denominator).undo_dots().to_expression()
        for position, momentum in (
            (0, (Q * n(idx) + Qbar * nb(idx)) / 2 + lam**2 * k(idx)),
            (1, Q * n(idx) / 2),
            (2, Qbar * nb(idx) / 2),
            (3, lam**2 * k(idx)),
        ):
            numerator = numerator.replace(P(position, idx), momentum)
            denominator = denominator.replace(P(position, idx), momentum)
        denominator = (
            kin.apply(
                TensorExpression(denominator).expand().simplify_metrics().to_dots()
            )
            .to_expression()
            .expand()
        )
        denominators.append(denominator)
        amplitude = (
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
            / denominator
        )
        leading = amplitude.series(
            lam, 0, -2 if len(outgoing) == 3 else 0
        ).to_expression()
        channels.append(simplify(project(leading.replace(lam, one))))
assert set(denominators) == {one, Q * kp * lam**2, Qbar * km * lam**2}
born = channels[0]
assert simplify(born * n(mink(D, ports[0]))) == zero
assert simplify(born * nb(mink(D, ports[0]))) == zero

# Keep the photon and both spinor ports open; divide out only the Born color delta.
color_identity = TensorExpression.g(
    Representation.cof(3), Representation.cof(3).dual()
)(ports[2], ports[1]).to_expression()
current = (born / color_identity).together()
generator = TensorExpression.t(8, 3)(ports[3], ports[2], ports[1]).to_expression()
gs = S("UFO::G")
# The generated convention gives a common minus relative to the reference's
# emission current. Its relative quark/antiquark sign and square agree.
for emitted, vector, product, sign in zip(channels[1:], (n, nb), (kp, km), (-1, 1)):
    expected = sign * gs * vector(mink(D, ports[3])) / product * current * generator
    assert simplify(emitted - expected).together() == zero
eikonal = ((channels[1] + channels[2]) / (current * generator * gs)).together()
assert (
    eikonal + n(mink(D, ports[3])) / kp - nb(mink(D, ports[3])) / km
).expand() == zero
assert simplify(eikonal * k(mink(D, ports[3]))) == zero
assert eikonal.derivative(Q) == zero and eikonal.derivative(Qbar) == zero
print(
    "PASS: both generated soft amplitudes factorize with open spinor/photon ports; Ward identity"
)

# Close color indices without summing or averaging the external spin states.
Nc, dA = S("soft::Nc", "soft::dA")
fund = Representation.cof(Nc)
color_born = TensorExpression.g(fund, fund.dual())(ports[2], ports[1])
color_real = TensorExpression.t(dA, Nc)(ports[3], ports[2], ports[1])
colors = []
for tensor in (color_born, color_real):
    norm = (
        TensorExpression(
            tensor.to_expression() * tensor.spenso_conjugate().to_expression()
        )
        .simplify_color()
        .to_expression()
    )
    colors.append(
        TensorExpression(norm.replace(dA, Nc**2 - 1))
        .to_cof_dimension_invariants()
        .to_expression()
    )
assert colors == [Nc, (Nc**2 - 1) / 2]
cf = (colors[1] / colors[0]).together()

# An arbitrary reference includes a nonzero r^2; no gauge term may survive.
rn, rnb, rk, r2 = S("soft::rn", "soft::rnb", "soft::rk", "soft::r2")
kin = (
    kin.with_scalar_product(ref, n, rn)
    .with_scalar_product(ref, nb, rnb)
    .with_scalar_product(ref, k, rk)
    .with_scalar_product(ref, ref, r2)
)
mu2 = S("soft::mu2")
polarization_results = {}
for label, reference in (
    ("Covariant", None),
    ("n", n),
    ("nbar", nb),
    ("Arbitrary reference", ref),
):
    density = model.particle("g").spin_sum(
        k, ports[3], mu2, dimension=D, reference=reference, covariant=reference is None
    )
    squared = simplify(eikonal * eikonal.replace(ports[3], mu2) * density).together()
    assert (squared - 4 / (kp * km)).together() == zero
    polarization_results[label] = squared
alpha = S("soft::alpha_s")
soft_function = (squared * cf * 4 * Symbol.PI * alpha / (32 * Symbol.PI**2)).together()
assert (soft_function - cf * alpha / (2 * Symbol.PI * kp * km)).together() == zero
eta = S("soft::eta")
assert (
    soft_function.replace(kp, eta * kp).replace(km, km / eta) - soft_function
).together() == zero
energy, z = S("soft::energy", "soft::z")
angular_soft = (
    soft_function.replace(kp, energy * (1 - z)).replace(km, energy * (1 + z)).together()
)
assert (angular_soft.replace(energy, 2 * energy) - angular_soft / 4).together() == zero
for nc in (2, 3, 5):
    for energy_value in (1, 2, 5):
        for angle in (-0.8, 0, 0.8):
            parameters = {Nc: nc, energy: energy_value, z: angle, alpha: 0.118}
            value = complex(angular_soft.evaluate(parameters)).real
            expected = (
                ((nc**2 - 1) / (2 * nc))
                * 0.118
                / (2 * 3.141592653589793 * energy_value**2 * (1 - angle**2))
            )
            assert abs(value - expected) < 1e-12
print(
    "PASS: four polarization sums, symbolic SU(N) soft function, basis rescaling and 27 numerical points"
)

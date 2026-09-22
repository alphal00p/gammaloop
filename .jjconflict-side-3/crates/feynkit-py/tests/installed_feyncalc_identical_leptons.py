"""Generated Bhabha/Moller interference and identical-electron event counting.

References:
https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-ElAel
https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElEl-ElEl

Use ordinary labeled amplitudes, retaining both diagrams and their ordering
signs. The gallery squares average initial spins and sum final spins; they do
not include the identical-final-state phase-space factor.
"""

from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
s, t, u, mass, charge = S(
    "identical_leptons::s",
    "identical_leptons::t",
    "identical_leptons::u",
    "UFO::Me",
    "UFO::ee",
)
kin = fk.Kinematics.mandelstam([P(i) for i in range(4)], [mass**2] * 4, [s, t, u])
ports = S(
    "identical_leptons::i0",
    "identical_leptons::i1",
    "identical_leptons::i2",
    "identical_leptons::i3",
)
a, b, c, inverse, wave, rep, index = S(
    "a_", "b_", "c_", "inverse_", "wave_", "rep_", "index_"
)
conjugate, wrapped = S("spenso::conj", "identical_leptons::adjoint_index")
pi = Expression.PI
s_cm = S("identical_leptons::s_cm", is_positive=True)
beta = S("identical_leptons::beta", is_positive=True)

# Equal incoming and outgoing masses give the same CM speed, 0<beta<1.
# The shared phase-space and flux factors therefore cancel this speed.
physical = fk.Kinematics.mandelstam(
    [P(i) for i in range(4)], [s_cm * (1 - beta**2) / 4] * 4, [s_cm, t, u]
)
flux = (
    physical.flux(P(0), P(1)).expand().replace((s_cm**2 * beta**2).sqrt(), s_cm * beta)
)
measure = (
    physical.two_body_phase_space(P(2), P(3))
    .expand()
    .replace((s_cm**2 * beta**2).sqrt(), s_cm * beta)
)
assert (flux - 2 * s_cm * beta).together() == E("0")
assert (measure - beta / (32 * pi**2)).together() == E("0")
assert (measure / flux - 1 / (64 * pi**2 * s_cm)).together() == E("0")

for label, pdgs in (
    ("Bhabha", [11, -11, 11, -11]),
    ("Moller", [11, 11, 11, 11]),
):
    generated = fk.Generator(model).generate(
        fk.Process.amplitude(pdgs[:2], pdgs[2:]),
        max_vertices=2,
        maximum_bridges=None,
        vertex_allow=["V_98"],
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 2
    operators, denominators, factors = [], [], []
    for diagram in generated.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
        for edge in diagram.external_edges:
            matches = list(
                diagram.projector_expression().match(
                    wave(edge.id, rep(4, index)), max_level=0
                )
            )
            assert len(matches) == 1
            numerator = numerator.replace(
                dict(matches[0])[index], ports[edge.external_index]
            )
        denominator = kin.apply(
            diagram.denominator_expression(dimension=4, in_lmb=True)
            .to_expression()
            .replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
        ).expand()
        denominators.append(denominator)
        factor = diagram.overall_factor_expression(evaluate=True)
        factors.append(factor)
        operator = TensorExpression(
            (
                numerator
                * factor
                * diagram.numerator_prefactor_expression()
                / denominator
            ).expand()
        )
        operators.append(operator.to_expression())
    assert set(denominators) == ({s, t} if label == "Bhabha" else {t, u})
    assert set(factors) == {E("1"), E("-1")}

    amplitude = sum(operators, E("0"))
    adjoints = []
    # Each tree diagram has two one-gamma currents with different pairings.
    # Preserve physical leg labels when adjoining the complete channel sum,
    # so common spin projectors close its direct and interference terms alike.
    for expression in [*operators, amplitude]:
        adjoint = (
            TensorExpression(expression)
            .dirac_adjoint(preserve_indices=True)
            .to_expression()
        )
        for real in (s, t, u, mass, charge):
            adjoint = adjoint.replace(conjugate(real), real)
        adjoints.append(TensorExpression(adjoint).wrap_indices(wrapped))
    combined_adjoint = adjoints.pop()
    assert (combined_adjoint - sum(adjoints, E("0"))).expand() == E("0")

    spins = E("1")
    for position, pdg in enumerate(pdgs):
        # Columns carry u or v; rows carry ubar or vbar. Average the two
        # incoming spins independently and sum both final spins.
        column = (position < 2) == (pdg > 0)
        left, right = (
            (ports[position], wrapped(ports[position]))
            if column
            else (wrapped(ports[position]), ports[position])
        )
        spins *= model.particle_by_pdg(pdg).spin_sum(
            P(position), left, right, average=position < 2
        )
    evaluated = []
    for raw in (
        amplitude * combined_adjoint,
        sum((x * y for x, y in zip(operators, adjoints, strict=True)), E("0")),
    ):
        scalar = (
            TensorExpression(raw * spins, cook_indices=CookSettings.indices())
            .expand()
            .simplify_gamma()
            .expand()
            .simplify_gamma()
            .simplify_metrics()
            .to_dots()
        )
        assert scalar.is_scalar
        evaluated.append(
            kin.apply(scalar.to_expression()).replace(u, 4 * mass**2 - s - t).together()
        )
    squared, diagonal = evaluated
    if label == "Bhabha":
        expected = (
            2
            * charge**4
            * (
                8 * mass**4 * (s**2 + s * t + t**2)
                - 4
                * mass**2
                * (s**3 + s**2 * (u - 2 * t) + s * t * (3 * u - 2 * t) + t**2 * (t + u))
                + s**4
                + s**2 * u**2
                + 2 * s * t * u**2
                + t**4
                + t**2 * u**2
            )
            / (s**2 * t**2)
        )
        expected_interference = (
            4 * charge**4 * (u**2 - 8 * mass**2 * u + 12 * mass**4) / (s * t)
        )
        expected_massless = (
            2
            * charge**4
            * ((s**2 + u**2) / t**2 + (t**2 + u**2) / s**2 + 2 * u**2 / (s * t))
        )
    else:
        expected = (
            2
            * charge**4
            * (
                -4
                * mass**2
                * (
                    s * (t**2 + 3 * t * u + u**2)
                    + t**3
                    - 2 * t**2 * u
                    - 2 * t * u**2
                    + u**3
                )
                + 8 * mass**4 * (t**2 + t * u + u**2)
                + s**2 * (t + u) ** 2
                + t**4
                + u**4
            )
            / (t**2 * u**2)
        )
        expected_interference = (
            4 * charge**4 * (s**2 - 8 * mass**2 * s + 12 * mass**4) / (t * u)
        )
        expected_massless = (
            2
            * charge**4
            * ((s**2 + u**2) / t**2 + (s**2 + t**2) / u**2 + 2 * s**2 / (t * u))
        )
        assert (squared - squared.replace(t, 4 * mass**2 - s - t)).together() == E("0")
    assert (squared - expected.replace(u, 4 * mass**2 - s - t)).together() == E("0"), (
        label
    )
    assert (
        squared - diagonal - expected_interference.replace(u, 4 * mass**2 - s - t)
    ).together() == E("0"), label
    massless = squared.replace(mass, E("0"))
    assert (massless - expected_massless.replace(u, -s - t)).together() == E("0"), label
    print(
        f"{label}: generated massive square, interference and massless limit passed",
        flush=True,
    )

    cos_theta, alpha, cutoff = S(
        "identical_leptons::cos_theta",
        "identical_leptons::alpha",
        "identical_leptons::cutoff",
    )
    angular = (
        massless.replace(t, -s * (1 - cos_theta) / 2)
        .replace(s, s_cm)
        .replace(charge**4, (4 * pi * alpha) ** 2)
    )
    labeled = (angular * measure / flux).together()
    kernel = (labeled * s_cm / alpha**2).together()
    expected_kernel = (3 + cos_theta**2) ** 2 / (
        4 * (1 - cos_theta) ** 2 if label == "Bhabha" else (1 - cos_theta**2) ** 2
    )
    assert (kernel - expected_kernel).together() == E("0")
    primitive = kernel.integrate(cos_theta).together()
    assert (primitive.derivative(cos_theta) - kernel).together() == E("0")
    # An angular cut 0<cutoff<1 excludes the forward/backward Coulomb poles.
    symmetric_integral = primitive.replace(cos_theta, cutoff) - primitive.replace(
        cos_theta, -cutoff
    )
    if label == "Bhabha":
        # Electron and positron are distinct: use the entire angular interval.
        full_sphere = 2 * pi * alpha**2 / s_cm * symmetric_integral
        expected_integral = (
            cutoff**3 / 6
            + 9 * cutoff / 2
            - 8 * cutoff.atanh()
            + 8 * cutoff / (1 - cutoff**2)
        )
        # Symbolica returns logarithms instead of atanh. Equality follows from
        # the same derivative and value at zero on the real interval |cutoff|<1.
        residual = (
            full_sphere * s_cm / (2 * pi * alpha**2) - expected_integral
        ).together()
        assert residual.derivative(cutoff).together() == E("0")
        assert residual.replace(cutoff, E("0")).together() == E("0")
        for value in ("1/5", "1/2", "4/5"):
            assert abs(complex(residual.replace(cutoff, E(value)).to_float(40))) < 1e-30
    else:
        assert (
            primitive - cos_theta - 8 * cos_theta / (1 - cos_theta**2)
        ).together() == E("0")
        # Either count one forward electron per event, or integrate labeled
        # electrons over the full symmetric interval with the explicit 1/2!.
        hemisphere = (
            2
            * pi
            * alpha**2
            / s_cm
            * (
                primitive.replace(cos_theta, cutoff)
                - primitive.replace(cos_theta, E("0"))
            )
        )
        full_sphere = E("1/2") * 2 * pi * alpha**2 / s_cm * symmetric_integral
        expected_events = (
            2 * pi * alpha**2 / s_cm * cutoff * (9 - cutoff**2) / (1 - cutoff**2)
        )
        assert (hemisphere - expected_events).together() == E("0")
        assert (full_sphere - expected_events).together() == E("0")
    print(
        f"{label}: angular-cut cross section and final-state event counting passed",
        flush=True,
    )

"""Generated massive photon-gluon currents in three crossed channels.

References:
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/GaGl-QQbar
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QGa-GlQ
https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QQbar-GaGl
"""

import numpy as np
from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import (
    CookSettings,
    GammaSimplifySettings,
    TensorExpression,
)

index_scope = S("spenso::index_scope")

model = hep.Model.standard_model()
P = S("gammalooprs::P")
mass, ee, gs = S("UFO::MB", "UFO::ee", "UFO::G")
s, t, u, virtuality = S(
    "photon_gluon::s", "photon_gluon::t", "photon_gluon::u", "photon_gluon::q2"
)
Nc, dA, cof, coad = S(
    "photon_gluon::Nc", "photon_gluon::dA", "spenso::cof", "spenso::coad"
)


def calculate(names, photon_position, gluon_position, fermion_ports):
    masses = [
        mass**2 if name in ("b", "b~") else virtuality if name == "a" else E("0")
        for name in names
    ]
    kin = hep.Kinematics.mandelstam([P(i) for i in range(4)], masses, [s, t, u])
    allowed = [
        v
        for v in model.vertex_rules
        if sorted(v.particles)
        in [
            sorted(["b", "b~", "a"]),
            sorted(["b", "b~", "g"]),
        ]
    ]
    generated = model.process(
        names[:2], names[2:], vertex_allow=allowed
    ).generate_diagrams(
        max_vertices=2, maximum_bridges=None, numerator_grouping=None, progress=None
    )
    assert len(generated.diagrams) == 2
    ports = S(
        "photon_gluon::i0", "photon_gluon::i1", "photon_gluon::i2", "photon_gluon::i3"
    )
    a, b, c, inv, index, left, right = S(
        "a_", "b_", "c_", "inv_", "index_", "left_", "right_"
    )
    conj, wrap, metric = S("spenso::conj", "photon_gluon::adjoint", "spenso::g")
    terms = []
    for diagram in generated.diagrams:
        numerator = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
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
        assert denominator in (s - mass**2, t - mass**2, u - mass**2), denominator
        terms.append(
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
            / denominator
        )
    operator = TensorExpression(sum(terms, E("0")))
    adjoint = (
        operator.dirac_adjoint()
        .simplify_gamma(GammaSimplifySettings(gamma0=True, evaluate_traces=False))
        .to_expression()
        .to_expression()
        .replace(conj(P(a, b)), P(a, b))
    )
    for real in (mass, ee, gs, s, t, u, virtuality):
        adjoint = adjoint.replace(conj(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(wrap).to_expression()
    color_projector = E("1")
    initial_colors = 1
    generic_initial_colors = E("1")
    for position, name in enumerate(names):
        particle = model.particle(name)
        if particle.color == 1:
            continue
        if position < 2:
            initial_colors *= abs(particle.color)
            generic_initial_colors *= dA if name == "g" else Nc
        matches = []
        for slot in operator.structure.slots:
            original = slot.to_expression()
            if original.replace(ports[position], E("0")) == original:
                continue
            closure = metric(
                slot.dual().to_expression(),
                original.replace(ports[position], index_scope(wrap, ports[position])),
            )
            match = next(
                closure.match(particle.color_sum(left, right), max_level=0), None
            )
            if match is not None:
                matches.append(dict(match))
        assert len(matches) == 1
        color_projector *= particle.color_sum(
            matches[0][left], matches[0][right], average=position < 2
        )
    generic = (
        (
            operator.to_expression()
            * adjoint
            * color_projector
            * initial_colors
            / generic_initial_colors
        )
        .replace(cof(3, index), cof(Nc, index))
        .replace(coad(8, index), coad(dA, index))
    )
    colored = (
        TensorExpression(generic, cook_indices=CookSettings.indices())
        .simplify_color()
        .to_expression()
        .to_expression()
        .replace(dA, Nc**2 - 1)
    )
    colored = TensorExpression(colored).to_cof_dimension_invariants()
    print("COLOR PASS", flush=True)
    # Adjoint reverses the single open fermion chain. The explicit endpoint
    # pairing distinguishes pair creation, annihilation and Compton scattering.
    spins = E("1")
    for position, left_position, right_position, wrap_left in fermion_ports:
        left_port, right_port = ports[left_position], ports[right_position]
        if wrap_left:
            left_port = index_scope(wrap, left_port)
        else:
            right_port = index_scope(wrap, right_port)
        spins *= model.particle(names[position]).spin_sum(
            P(position), left_port, right_port, average=position < 2
        )
    traced = (
        TensorExpression(
            colored.to_expression() * spins,
            cook_indices=CookSettings.indices(),
        )
        .simplify_gamma()
        .to_expression()
    )
    reduced = kin.apply(traced.contract().to_expression().to_dots())
    results = {}
    modes = (
        "covariant",
        "photon reference",
        "quark reference",
        "gluon Ward",
        "photon Ward",
    )
    quark_position = next(i for i, name in enumerate(names) if name == "b")
    for mode in modes:
        photon = model.particle("a").spin_sum(
            P(photon_position),
            ports[photon_position],
            index_scope(wrap, ports[photon_position]),
            covariant=True,
        )
        reference = (
            P(photon_position)
            if mode == "photon reference"
            else P(quark_position)
            if mode == "quark reference"
            else None
        )
        gluon = model.particle("g").spin_sum(
            P(gluon_position),
            ports[gluon_position],
            index_scope(wrap, ports[gluon_position]),
            reference=reference,
            covariant=reference is None,
            average=gluon_position < 2,
        )
        if mode == "gluon Ward":
            gluon = P(gluon_position, S("spenso::mink")(4, ports[gluon_position])) * P(
                gluon_position,
                S("spenso::mink")(4, index_scope(wrap, ports[gluon_position])),
            )
        if mode == "photon Ward":
            photon = P(
                photon_position, S("spenso::mink")(4, ports[photon_position])
            ) * P(
                photon_position,
                S("spenso::mink")(4, index_scope(wrap, ports[photon_position])),
            )
        # The shared contractor selects the connected tensor factors while
        # retaining scalar coefficients; no coefficient-splitting API is needed.
        contracted = (
            TensorExpression(
                reduced.to_expression() * photon * gluon,
                cook_indices=CookSettings.indices(),
            )
            .contract()
            .to_expression()
            .to_dots()
        )
        assert contracted.is_scalar
        results[mode] = (
            kin.apply(contracted.to_expression())
            .replace(u, 2 * mass**2 + virtuality - s - t)
            .together()
        )
        print("PASS", names, mode, flush=True)
    for mode in [name for name in modes if name.endswith("Ward")]:
        assert results[mode] == E("0")
    for mode in [name for name in modes if name.endswith("reference")]:
        assert (results[mode] - results["covariant"]).together() == E("0")
    print("ALL GAUGES PASS", flush=True)

    return generated, results, kin


# Reference polynomial shared by the crossed channels, independently stated
# in invariant variables. No crossing rule is inserted into the amplitudes.
def reference_polynomial(x, y):
    return (
        -(mass**4)
        * (
            2 * virtuality**2
            - 2 * virtuality * (x + y)
            + 3 * x**2
            + 14 * x * y
            + 3 * y**2
        )
        + mass**2
        * (
            2 * virtuality**2 * (x + y)
            - 8 * virtuality * x * y
            + x**3
            + 7 * x**2 * y
            + 7 * x * y**2
            + y**3
        )
        + 6 * mass**8
        - x * y * (2 * virtuality**2 - 2 * virtuality * (x + y) + x**2 + y**2)
    ) / ((x - mass**2) ** 2 * (y - mass**2) ** 2)


channels = {}
for channel, names, photon, gluon, spin_ports, prefactor, pair in (
    (
        "fusion",
        ["a", "g", "b", "b~"],
        0,
        1,
        [(2, 3, 2, True), (3, 3, 2, False)],
        -2 * ee**2 * gs**2 / 9,
        (t, u),
    ),
    (
        "compton",
        ["b", "a", "g", "b"],
        1,
        2,
        [(0, 0, 3, False), (3, 0, 3, True)],
        2 * ee**2 * gs**2 * (Nc**2 - 1) / (9 * Nc),
        (s, t),
    ),
    (
        "annihilation",
        ["b", "b~", "a", "g"],
        2,
        3,
        [(0, 0, 1, False), (1, 0, 1, True)],
        -(ee**2) * gs**2 * (Nc**2 - 1) / (9 * Nc**2),
        (t, u),
    ),
):
    generated, results, kin = calculate(names, photon, gluon, spin_ports)
    target = (
        (prefactor * reference_polynomial(*pair))
        .replace(u, 2 * mass**2 + virtuality - s - t)
        .together()
    )
    delta = (results["covariant"] - target).together()
    print("REFERENCE PASS", channel, flush=True)
    assert delta == E("0"), channel
    x, y = pair
    massless_target = (
        -prefactor
        * (2 * virtuality**2 - 2 * virtuality * (x + y) + x**2 + y**2)
        / (x * y)
    )
    assert (
        results["covariant"].replace(mass, 0)
        - massless_target.replace(u, virtuality - s - t)
    ).together() == E("0")
    channels[channel] = (generated, results, kin)
print(
    "PASS: all three massive currents, physical gluon gauges, Ward identities and massless limits",
    flush=True,
)

fusion = channels["fusion"][1]["covariant"]
compton = channels["compton"][1]["covariant"]
annihilation = channels["annihilation"][1]["covariant"]
# Crossing changes both the fermion-trace sign and the incoming color average.
crossed_fusion = fusion.replace(s, 2 * mass**2 + virtuality - s - t)
assert (compton + (Nc**2 - 1) / Nc * crossed_fusion).together() == E("0")
assert (annihilation - (Nc**2 - 1) / (2 * Nc**2) * fusion).together() == E("0")
assert (fusion - fusion.replace(t, 2 * mass**2 + virtuality - s - t)).together() == E(
    "0"
)
# Check the reference spacelike continuation q^2=-Q^2 explicitly at zero mass.
Q2 = S("photon_gluon::Q2", is_positive=True)
for channel, x, y, prefactor in (
    ("fusion", t, u, 2 * ee**2 * gs**2 / 9),
    ("compton", s, t, -16 * ee**2 * gs**2 / 27),
):
    reference = prefactor * (x / y + y / x + 2 * Q2 * (x + y + Q2) / (x * y))
    result = (
        channels[channel][1]["covariant"]
        .replace(mass, 0)
        .replace(virtuality, -Q2)
        .replace(Nc, 3)
    )
    assert (result - reference.replace(u, -Q2 - s - t)).together() == E("0")

# Only for a real photon do we form ordinary two-body cross sections.
# Incoming photons have two physical states, so include their missing 1/2;
# the off-shell current convention used above intentionally has no such average.
z, rho, beta = S("photon_gluon::z", "photon_gluon::rho", "photon_gluon::beta")
alpha, alpha_s = S("photon_gluon::alpha", "photon_gluon::alpha_s")
angular = {}
for channel, (_, sums, kin) in channels.items():
    phase = (
        (kin.two_body_phase_space(P(2), P(3)) / kin.flux(P(0), P(1)))
        .replace(virtuality, 0)
        .together()
    )
    # Squaring the measure ratio checks its invariant normalization without
    # imposing a branch identity outside the physical region.
    velocity_squared = 1 - 4 * mass**2 / s
    factor_squared = (
        velocity_squared
        if channel == "fusion"
        else 1 / velocity_squared
        if channel == "annihilation"
        else E("1")
    )
    assert (phase**2 - factor_squared / (64 * Symbol.PI**2 * s) ** 2).together() == E(
        "0"
    )
    angle_t = (
        mass**2 - s * (1 - beta * z) / 2
        if channel != "compton"
        else mass**2 - (s**2 - mass**4) / (2 * s) + (s - mass**2) ** 2 * z / (2 * s)
    )
    # Normalize out alpha*alpha_s/s. Final particles are distinct in all channels.
    density = (
        sums["covariant"].replace(virtuality, 0).replace(t, angle_t)
        * phase
        * 2
        * Symbol.PI
        * s
        / (alpha * alpha_s)
    )
    density = density.replace(ee**2, 4 * Symbol.PI * alpha).replace(
        gs**2, 4 * Symbol.PI * alpha_s
    )
    if channel != "annihilation":
        density /= 2
    density = density.replace(mass, (rho * s).sqrt()).together()
    angular[channel] = density
    for nc in (2, 3, 5):
        for r in (0.0, 0.02, 0.1):
            for cosine in (-0.7, 0.0, 0.7):
                values = {s: 100.0, Nc: nc, rho: r, beta: np.sqrt(1 - 4 * r), z: cosine}
                exact = complex(density.evaluate(values))
                scaled = complex(density.evaluate({**values, s: 400.0}))
                assert abs(exact.imag) < 1e-12 and exact.real > 0
                assert abs(exact - scaled) < 1e-10

# Off-shell currents at explicit center-of-mass invariants, with incoming
# photons spacelike and the outgoing photon timelike.
for channel, (_, sums, _) in channels.items():
    for r in (0.0, 0.02, 0.1):
        for v in (0.05, 0.2):
            v = v if channel == "annihilation" else -v
            for cosine in (-0.7, 0.0, 0.7):
                kallen = 1 + r * r + v * v - 2 * r - 2 * v - 2 * r * v
                angle_t = (
                    r - (1 - v) * (1 - np.sqrt(1 - 4 * r) * cosine) / 2
                    if channel != "compton"
                    else r
                    - (1 + r - v) * (1 - r) / 2
                    + (1 - r) * np.sqrt(kallen) * cosine / 2
                )
                point = {
                    s: 1,
                    t: angle_t,
                    mass: np.sqrt(r),
                    virtuality: v,
                    Nc: 3,
                    ee: 1,
                    gs: 1,
                }
                values = [
                    complex(sums[mode].evaluate(point))
                    for mode in ("covariant", "photon reference", "quark reference")
                ]
                assert max(abs(value.imag) for value in values) < 1e-10
                assert max(abs(value - values[0]) for value in values) < 1e-8
                assert all(np.isfinite(value.real) for value in values)
print(
    "PASS: exact crossing, spacelike continuation, 54 off-shell points and 81 on-shell angular/scaling checks",
    flush=True,
)

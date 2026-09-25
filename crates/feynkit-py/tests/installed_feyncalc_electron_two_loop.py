"""Generated massless QED electron self-energy in Feynman gauge (xi=1).

Primary reference:
https://feyncalc.github.io/FeynCalcExamples/QED/TwoLoops/Renormalization-LeAle-Massless
The equal-mass vacuum master poles and one-loop counterterm insertion sum are
explicit analytic inputs. Actual generated diagrams, UV expansion, tensor
projection, topology mappings and bounded native IBP produce the bare poles.
"""

from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
d, s, mass_squared, uv_mass, eps, log_mass, log_4pi, flavors, coupling = S(
    "electron_2l::D",
    "electron_2l::s",
    "electron_2l::M",
    "electron_2l::mUV",
    "electron_2l::eps",
    "electron_2l::Lm",
    "electron_2l::L4pi",
    "electron_2l::Nf",
    "electron_2l::a4",
)
k, p, mass, charge = S("gammalooprs::K", "gammalooprs::P", "UFO::Me", "UFO::ee")
mink, bis, gamma, wave = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "electron_2l::wave_"
)
index, dimension, mu = S(
    "electron_2l::index_", "electron_2l::dimension_", "electron_2l::mu"
)
ordering, closed_loop, value = S(
    "feynkit_generator_factor::ExternalFermionOrderingSign",
    "feynkit_generator_factor::InternalFermionLoopSign",
    "electron_2l::value_",
)
denominator, edge_pattern, momentum_pattern, mass_pattern, quadratic_pattern = S(
    "gammalooprs::denom",
    "electron_2l::edge_",
    "electron_2l::momentum_",
    "electron_2l::mass_",
    "electron_2l::quadratic_",
)
coordinates = list(S("electron_2l::d0", "electron_2l::d1", "electron_2l::d2"))
tadpole_squared, sunset, integral = S(
    "electron_2l::T2", "electron_2l::V", "electron_2l::I"
)
zero, one, imaginary = E("0"), E("1"), E("1i")
kinematics = hep.Kinematics(d, momenta=[k(0), k(1), p(0)]).with_scalar_product(
    p(0), p(0), s
)
vacuum_kinematics = hep.Kinematics(d, momenta=[k(0), k(1)])
family = hep.IntegralFamily(
    [k(0), k(1)],
    [],
    [
        vacuum_kinematics.scalar_product(q, q) - mass_squared
        for q in (k(0), k(1), k(0) - k(1))
    ],
    kinematics=vacuum_kinematics,
)
reducer = (
    hep.TensorReducer(d)
    .with_integrated_vector(k(0, mink(d)))
    .with_integrated_vector(k(1, mink(d)))
)

# Pin the 1PI convention independently at one loop, before any UV expansion.
# Vertices (-ie)^2, electron i*slash(k), photon -i*g give -e^2 gamma_mu slash(k) gamma^mu.
# Native amplitudes additionally order external fields by leg ID: [0,1] rather
# than the conventional kernel order [antifermion,fermion]=[1,0].
one_loop = (
    hep.Process(model, [11], [11])
    .with_loop_count(1, 1)
    .generate_diagrams(
        max_vertices=2,
        maximum_bridges=0,
        vertex_allow=["V_98"],
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    .diagrams[0]
)
assert one_loop.overall_factor_expression().matches(ordering(-1))
one_ports = {
    edge.external_index: dict(
        next(
            one_loop.projector_expression().match(
                wave(edge.id, bis(4, index)), max_level=0
            )
        )
    )[index]
    for edge in one_loop.external_edges
}
one_numerator = model.expand_couplings(
    one_loop.numerator_expression(in_lmb=True).to_expression()
).replace(mass, zero)
one_numerator = one_numerator.replace(mink(4, index), mink(d, index))
one_numerator *= (
    gamma(bis(4, one_ports[0]), bis(4, one_ports[1]), mink(d, mu))
    * p(0, mink(d, mu))
    / (4 * s)
)
one_trace = kinematics.apply(
    TensorExpression(one_numerator.expand())
    .simplify_gamma()
    .expand()
    .to_dots()
    .to_expression()
)
assert (
    one_trace + charge**2 * (2 - d) * kinematics.scalar_product(k(0), p(0)) / s
).together() == zero
one_kinematics = hep.Kinematics(d, momenta=[k(0), p(0)]).with_scalar_product(
    p(0), p(0), s
)
one_family = hep.IntegralFamily(
    [k(0)],
    [p(0)],
    [one_kinematics.scalar_product(q, q) for q in (k(0), k(0) - p(0))],
    kinematics=one_kinematics,
)
reflection = one_family.mapping_to(one_family, [p(0) - k(0)])
assert reflection is not None
assert (
    one_trace + reflection.apply(one_trace) + charge**2 * (2 - d)
).together() == zero
assert abs(complex(oneloop.B0(-1.0, 0.0, 0.0, 1.0)[1]) - 1) < 1e-12
# Thus the normalized one-loop kernel has +i*a4*slash(p)/eps and deltaZpsi=-1/eps.
zpsi_one = -one / eps

result = (
    hep.Process(model, [11], [11])
    .with_loop_count(2, 2)
    .generate_diagrams(
        max_vertices=4,
        maximum_bridges=0,
        vertex_allow=["V_98"],
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
)
assert len(result.diagrams) == 3
targets, diagram_integrals, flavor_weights = set(), [], []
for diagram in result.diagrams:
    raw_factor = diagram.overall_factor_expression()
    kernel_factor = raw_factor.replace(ordering(value), one)
    removed_ordering = (raw_factor / kernel_factor).replace(ordering(value), value)
    assert removed_ordering == -one
    # Evaluate every other native factor unchanged. In particular, a closed
    # fermion loop keeps its minus sign when external Wick ordering is removed.
    factor = diagram.overall_factor_expression(evaluate=True) / removed_ordering
    has_closed_loop = bool(raw_factor.matches(closed_loop(-1)))
    assert bool(kernel_factor.matches(closed_loop(-1))) == has_closed_loop
    assert factor == (-one if has_closed_loop else one)
    flavor_weights.append(flavors if has_closed_loop else one)
    numerator = model.expand_couplings(
        diagram.numerator_expression().to_expression()
    ).replace(mass, zero)
    uv = diagram.momentum_basis().route_expression(
        diagram.uv_expansion(uv_mass, numerator=numerator).to_expression()
    )
    uv = (
        uv.replace(mink(dimension, index), mink(d, index))
        .replace(mink(dimension), mink(d))
        .replace(uv_mass**2, mass_squared)
    )
    ports = {
        edge.external_index: dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, bis(4, index)), max_level=0
                )
            )
        )[index]
        for edge in diagram.external_edges
    }
    uv *= (
        gamma(bis(4, ports[0]), bis(4, ports[1]), mink(d, mu))
        * p(0, mink(d, mu))
        / (4 * s)
    )
    pattern = denominator(
        edge_pattern, momentum_pattern, mass_pattern, quadratic_pattern
    )
    quadratics = {dict(match)[quadratic_pattern] for match in uv.match(pattern)}
    source_family = hep.IntegralFamily(
        [k(0), k(1)], [], sorted(quadratics, key=str), kinematics=vacuum_kinematics
    )
    mapping = source_family.find_mapping(family)
    assert mapping is not None
    # Keep inverse denominators opaque during the polynomial tensor projection.
    for match in list(uv.match(pattern)):
        values = dict(match)
        formal = source_family.rewrite_numerator(values[quadratic_pattern], coordinates)
        assert formal in coordinates
        uv = uv.replace(
            denominator(
                values[edge_pattern],
                values[momentum_pattern],
                values[mass_pattern],
                values[quadratic_pattern],
            ),
            formal,
        )
    traced = (
        TensorExpression(uv.expand())
        .simplify_gamma()
        .expand()
        .to_dots()
        .to_expression()
    )
    scalar = source_family.rewrite_numerator(
        kinematics.apply(reducer.reduce(traced)), coordinates
    )
    scalar *= factor * diagram.numerator_prefactor_expression()
    terms = []
    for monomial, coefficient in scalar.expand().coefficient_list(*coordinates):
        powers = [
            -int((monomial.derivative(den) * den / monomial).together())
            for den in coordinates
        ]
        assert (
            monomial
            == coordinates[0] ** -powers[0]
            * coordinates[1] ** -powers[1]
            * coordinates[2] ** -powers[2]
        )
        powers = mapping.map_powers(powers)
        targets.add(tuple(powers))
        terms.append((powers, coefficient))
    diagram_integrals.append(terms)
assert flavor_weights.count(flavors) == 1
assert len(targets) == 22
solution = hep.IBPFamily(family, name="electron_two_loop").reduce_laporta(
    [list(target) for target in sorted(targets)], max_depth=2
)
assert solution.stats["rows"] > 0
assert {tuple(powers) for powers in solution.residuals} == {
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
    (1, 1, 1),
}
# Verified unit-Jacobian maps relate the three disconnected tadpole products.
for images in ([k(1), k(0)], [k(0), k(0) - k(1)]):
    assert family.mapping_to(family, images) is not None
basis = {
    (0, 1, 1): tadpole_squared,
    (1, 0, 1): tadpole_squared,
    (1, 1, 0): tadpole_squared,
    (1, 1, 1): sunset,
}
basis_rules = [
    Replacement(integral(*powers), master) for powers, master in basis.items()
]
# Explicit analytic inputs, not values inferred or certified by the IBP solver.
master_poles = [
    Replacement(
        tadpole_squared, mass_squared**2 * (1 / eps**2 + 2 * (1 - log_mass) / eps)
    ),
    Replacement(
        sunset, mass_squared * (E("3/2") / eps**2 + (E("9/2") - 3 * log_mass) / eps)
    ),
]
bare_uv = zero
for terms, flavor_weight in zip(diagram_integrals, flavor_weights, strict=True):
    reduced = sum(
        (
            coefficient * solution.reduce(powers, integral=integral)
            for powers, coefficient in terms
        ),
        zero,
    )
    reduced = reduced.replace_multiple(basis_rules).together()
    poles = reduced.replace_multiple(master_poles).replace(d, 4 - 2 * eps)
    # Each loop measure is i*(4pi)^(eps-2); divide by i*e^4/(16pi^2)^2.
    poles = (
        (-poles * (1 + 2 * eps * log_4pi) / (imaginary * charge**4))
        .series(eps, 0, -1)
        .to_expression()
    )
    bare_uv += flavor_weight * poles
expected_bare = (
    one / (2 * eps**2) + (log_4pi - log_mass - E("17/12") - E("7/3") * flavors) / eps
)
assert (bare_uv - expected_bare).together() == zero

# Reference one-loop CT insertion sum, in the same i*a4^2*slash(p) units.
# Includes the auxiliary photon-mass CT: deltaZAm=-2*Nf/eps; the other
# one-loop inputs are deltaZA=deltaZxi=-4*Nf/(3eps), deltaZe=2*Nf/(3eps),
# deltaZpsi=-1/eps. Automatic forest/CT generation is not claimed here.
counterterm_uv = -one / eps**2 + (4 * flavors / 3 + E("2/3") - log_4pi + log_mass) / eps
zpsi_two = -(bare_uv + counterterm_uv).expand()
assert zpsi_two.derivative(log_mass).expand() == zero
assert zpsi_two.derivative(log_4pi).expand() == zero
assert (
    zpsi_two - one / (2 * eps**2) - (4 * flavors + 3) / (4 * eps)
).together() == zero
zpsi = 1 + coupling * zpsi_one + coupling**2 * zpsi_two
assert (
    zpsi
    - 1
    + coupling / eps
    - coupling**2 * (one / (2 * eps**2) + (4 * flavors + 3) / (4 * eps))
).together() == zero
print(
    "Electron two-loop xi=1: generated bare UV, preserved closed-loop sign, 22-target IBP and reference-CT Zpsi passed",
    solution.stats,
)

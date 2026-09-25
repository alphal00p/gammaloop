"""Generated two-loop massless QED photon UV and calculated one-loop insertions.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/TwoLoops/Renormalization-GaGa
Use the reference's selective auxiliary-mass rearrangement, with symbolic gauge
parameter and tr(1)=4. Equal-mass vacuum poles, the tadpole finite part, and the
one-loop deltaZpsi (with deltaZ1=deltaZpsi) are explicit inputs. The four local
CT insertions and their integrated sum are calculated from a generated bubble.
No automatic counterterm diagrams, subtraction forest or finite two-loop
self-energy is claimed; bounded IBP residuals are not certified masters.
"""

import json
from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
d, s, mass_squared, uv_mass, eps, log_mass, log_4pi, flavors = S(
    "photon_2l::D",
    "photon_2l::s",
    "photon_2l::M",
    "photon_2l::mUV",
    "photon_2l::eps",
    "photon_2l::Lm",
    "photon_2l::L4pi",
    "photon_2l::Nf",
)
k, p, mass, charge = S("gammalooprs::K", "gammalooprs::P", "UFO::Me", "UFO::ee")
mink, wave = S("spenso::mink", "photon_2l::wave_")
index, dimension, mu = S("photon_2l::index_", "photon_2l::dimension_", "photon_2l::mu")
closed_loop, arguments = S(
    "feynkit_generator_factor::InternalFermionLoopSign",
    "photon_2l::arguments___",
)
denominator, edge_pattern, momentum_pattern, mass_pattern, quadratic_pattern = S(
    "gammalooprs::denom",
    "photon_2l::edge_",
    "photon_2l::momentum_",
    "photon_2l::mass_",
    "photon_2l::quadratic_",
)
coordinates = list(S("photon_2l::d0", "photon_2l::d1", "photon_2l::d2"))
tadpole_squared, sunset, integral = S("photon_2l::T2", "photon_2l::V", "photon_2l::I")
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


xi, nu = S("photon_2l::xi", "photon_2l::nu")
metric = S("spenso::g")
specification = json.loads(model.to_json())
for particle in specification["particles"]:
    if abs(particle["pdg_code"]) == 11:
        particle["mass"] = "ZERO"
for propagator in specification["propagators"]:
    if propagator["particle"] == "a":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))-(1-photon_2l::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(specification))
gmunu = metric(mink(d, mu), mink(d, nu))
ppmunu = p(0, mink(d, mu)) * p(0, mink(d, nu))
Qg, Qpp = S("photon_2l::metric_basis", "photon_2l::momentum_basis")
Q, dot = S("gammalooprs::Q", "spenso::dot")
result = (
    hep.Process(model, [22], [22])
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
targets, diagram_integrals = set(), []
for diagram in result.diagrams:
    raw_factor = diagram.overall_factor_expression()
    assert raw_factor.matches(closed_loop(-one))
    factor = diagram.overall_factor_expression(evaluate=True)
    numerator = model.expand_couplings(
        diagram.numerator_expression().to_expression()
    ).replace(mass, zero)
    photon = next(edge for edge in diagram.internal_edges if edge.particle_pdg == 22)
    # Selective mass replacement acts on the two propagator terms separately.
    # Clearing q² across the Feynman term would instead introduce an extra
    # M/(q²-M)² contribution and change this infrared-rearrangement prescription.
    feynman = numerator.replace(xi, one)
    longitudinal = (
        ((numerator - feynman) * dot(Q(photon.id, mink(4)), Q(photon.id, mink(4))))
        .together()
        .expand()
    )
    uv = diagram.uv_expansion(uv_mass, numerator=feynman).to_expression()
    uv += diagram.uv_expansion(
        uv_mass, numerator=longitudinal, edge_powers={photon.id: 2}
    ).to_expression()
    uv = diagram.momentum_basis().route_expression(uv)
    uv = uv.replace(mink(dimension, index), mink(d, index)).replace(
        mink(dimension), mink(d)
    )
    ports = {
        edge.external_index: dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, mink(4, index)), max_level=0
                )
            )
        )[index]
        for edge in diagram.external_edges
    }
    uv = uv.replace(mink(d, ports[0]), mink(d, mu)).replace(
        mink(d, ports[1]), mink(d, nu)
    )
    pattern = denominator(
        edge_pattern, momentum_pattern, mass_pattern, quadratic_pattern
    )
    quadratics = {
        dict(match)[quadratic_pattern].replace(uv_mass**2, mass_squared)
        for match in uv.match(pattern)
    }
    source_family = hep.IntegralFamily(
        [k(0), k(1)], [], sorted(quadratics, key=str), kinematics=vacuum_kinematics
    )
    mapping = source_family.find_mapping(family)
    assert mapping is not None
    # Freeze the massive denominators before removing explicit mUV terms.
    # Selecting mUV^0 discards the compensation terms in the shared full UV
    # expansion while retaining M in the integral family.  This reproduces the
    # reference's selective auxiliary-mass replacement and momentum Taylor series;
    # it does not change the default full-UV prescription used by other notebooks.
    # Keep inverse denominators opaque during the polynomial tensor projection.
    for match in list(uv.match(pattern)):
        values = dict(match)
        formal = source_family.rewrite_numerator(
            values[quadratic_pattern].replace(uv_mass**2, mass_squared), coordinates
        )
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
    assert not uv.matches(pattern)
    for monomial, _ in uv.expand().coefficient_list(uv_mass):
        power = int((monomial.derivative(uv_mass) * uv_mass / monomial).together())
        assert power >= 0 and monomial == uv_mass**power
    uv = uv.replace(uv_mass, zero)
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
    # Scalar basis labels avoid differentiating symbolic-D tensor slots during
    # the epsilon expansion.  The complete open tensor is checked before this.
    cg = scalar.expand().coefficient(gmunu)
    cp = scalar.expand().coefficient(ppmunu)
    assert (scalar - cg * gmunu - cp * ppmunu).expand() == zero
    scalar = cg * Qg + cp * Qpp
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
        assert all(coefficient.derivative(x).expand() == zero for x in coordinates), (
            coefficient
        )
        assert not coefficient.matches(k(arguments))
        assert not coefficient.matches(p(arguments))
        assert not coefficient.matches(Q(arguments))
        powers = mapping.map_powers(powers)
        targets.add(tuple(powers))
        terms.append((powers, coefficient))
    diagram_integrals.append(terms)
solution = hep.IBPFamily(family, name="photon_two_loop").reduce_laporta(
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
diagram_poles = []
for terms in diagram_integrals:
    reduced = sum(
        (
            coefficient * solution.reduce(powers, integral=integral)
            for powers, coefficient in terms
        ),
        zero,
    )
    reduced = reduced.replace_multiple(basis_rules).together()
    # Pole-only analytic masters suffice only if their exact coefficients have
    # no spurious pole at D=4.  The two independent master symbols stay intact.
    assert reduced.replace(d, 4 - 2 * eps).series(eps, 0, -1).to_expression() == zero, (
        "Master coefficient singular at D=4"
    )
    poles = reduced.replace_multiple(master_poles).replace(d, 4 - 2 * eps)
    # Each loop measure is i*(4pi)^(eps-2); divide by i*e^4/(16pi^2)^2.
    poles = (
        (-poles * (1 + 2 * eps * log_4pi) / (imaginary * charge**4))
        .series(eps, 0, -1)
        .to_expression()
    )
    poles = poles.expand()
    diagram_poles.append(poles)
    bare_uv += flavors * poles

bare_uv = bare_uv.expand()

transverse = s * Qg - Qpp
expected_bare = flavors * (
    mass_squared * xi * Qg * (-2 / eps**2 + (4 * log_mass - 4 * log_4pi + 3) / eps)
    - 2 * (2 * xi + 3) * transverse / (3 * eps)
)
assert (bare_uv - expected_bare).together() == zero

ct_coordinate, ct_integral = S("photon_2l::ctd", "photon_2l::ctI")
# Derive the four one-loop counterterm insertions on one generated bubble.
# QED.mod gives vertexCT/tree = sqrt(ZA)*Ze*Zpsi-1 and a kinetic insertion
# i*(Zpsi-1)*slash(q).  The one-loop Ward relation gives deltaZ1=deltaZpsi;
# hence two vertex factors +deltaZpsi and two kinetic factors -deltaZpsi*q_i².
# Keeping the latter's squared denominator through IR rearrangement is essential.
# Physical electron mass is zero, so its mass counterterm contributes nothing.
# These explicit local insertions are not automatic CT/forest generation.
one_loop_generated = (
    hep.Process(model, [22], [22])
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
)
assert len(one_loop_generated.diagrams) == 1
one_loop_diagram = one_loop_generated.diagrams[0]
factor = (
    one_loop_diagram.overall_factor_expression(evaluate=True)
    * one_loop_diagram.numerator_prefactor_expression()
)
assert factor == -one  # Preserve the native closed-fermion-loop sign.
fermions = [edge for edge in one_loop_diagram.edges if not edge.is_external]
assert len(fermions) == 2 and all(abs(edge.particle_pdg) == 11 for edge in fermions)
ports = {
    edge.external_index: dict(
        next(
            one_loop_diagram.projector_expression().match(
                wave(edge.id, mink(4, index)), max_level=0
            )
        )
    )[index]
    for edge in one_loop_diagram.external_edges
}
numerator = model.expand_couplings(
    one_loop_diagram.numerator_expression().to_expression()
).replace(mass, zero)
ct_kinematics = hep.Kinematics(d, momenta=[k(0), p(0)]).with_scalar_product(
    p(0), p(0), s
)
ct_vacuum = hep.Kinematics(d, momenta=[k(0)])
ct_family = hep.IntegralFamily(
    [k(0)],
    [],
    [ct_vacuum.scalar_product(k(0), k(0)) - mass_squared],
    kinematics=ct_vacuum,
)
ct_reducer = hep.TensorReducer(d).with_integrated_vector(k(0, mink(d)))
# Supplied one-loop field counterterm; the bubble and insertion sum are computed.
zpsi_one = -xi / eps
insertions = [
    ("bare", numerator, {}, one),
    ("vertex_0", numerator, {}, zpsi_one),
    ("vertex_1", numerator, {}, zpsi_one),
]
for number, edge in enumerate(fermions):
    qi = Q(edge.id, mink(4))
    insertions.append(
        (f"kinetic_{number}", numerator * dot(qi, qi), {edge.id: 2}, -zpsi_one)
    )
ct_parts, ct_targets = ({}, set())
for label, inserted, powers, weight in insertions:
    uv = one_loop_diagram.momentum_basis().route_expression(
        one_loop_diagram.uv_expansion(
            uv_mass, numerator=inserted, edge_powers=powers
        ).to_expression()
    )
    uv = uv.replace(mink(dimension, index), mink(d, index)).replace(
        mink(dimension), mink(d)
    )
    uv = uv.replace(mink(d, ports[0]), mink(d, mu)).replace(
        mink(d, ports[1]), mink(d, nu)
    )
    pattern = denominator(
        edge_pattern, momentum_pattern, mass_pattern, quadratic_pattern
    )
    for match in list(uv.match(pattern)):
        values = dict(match)
        coordinate = ct_family.rewrite_numerator(
            values[quadratic_pattern].replace(uv_mass**2, mass_squared), [ct_coordinate]
        )
        assert coordinate == ct_coordinate
        uv = uv.replace(
            denominator(
                values[edge_pattern],
                values[momentum_pattern],
                values[mass_pattern],
                values[quadratic_pattern],
            ),
            coordinate,
        )
    assert not uv.matches(pattern)
    for monomial, _ in uv.expand().coefficient_list(uv_mass):
        power = int((monomial.derivative(uv_mass) * uv_mass / monomial).together())
        assert power >= 0 and monomial == uv_mass**power
    # As above, formal coordinates already retain the massive denominators.
    selective = uv.replace(uv_mass, zero)
    trace = (
        TensorExpression(selective.expand())
        .simplify_gamma()
        .expand()
        .to_dots()
        .to_expression()
    )
    scalar = ct_family.rewrite_numerator(
        ct_kinematics.apply(ct_reducer.reduce(trace)), [ct_coordinate]
    )
    tensor = (factor * scalar / charge**2).together().expand()
    scalar = tensor.replace(gmunu, Qg).replace(ppmunu, Qpp).expand()
    assert (
        scalar - scalar.coefficient(Qg) * Qg - scalar.coefficient(Qpp) * Qpp
    ).expand() == zero
    terms = []
    for monomial, coefficient in scalar.coefficient_list(ct_coordinate):
        power = -int(
            (monomial.derivative(ct_coordinate) * ct_coordinate / monomial).together()
        )
        assert monomial == ct_coordinate ** (-power)
        assert coefficient.derivative(ct_coordinate).expand() == zero
        assert not coefficient.matches(k(arguments))
        assert not coefficient.matches(p(arguments))
        assert not coefficient.matches(Q(arguments))
        ct_targets.add((power,))
        terms.append(([power], coefficient))
    ct_parts[label] = terms
ct_solution = hep.IBPFamily(ct_family, name="photon2_ct").reduce_laporta(
    [list(target) for target in sorted(ct_targets)], max_depth=2
)
assert ct_solution.residuals == [[1]]
# The analytic A0 finite term is an explicit input: multiplication by deltaZpsi
# promotes it to a single pole.  Retain the D dependence and loop-measure term.
# ct_sum is divided by i*a4²*Nf; bare_uv already includes the flavor multiplicity.
tadpole = mass_squared * (1 / eps + 1 - log_mass)
ct_reduced, ct_integrated, ct_poles = ({}, {}, {})
for label, _, _, weight in insertions:
    ct_reduced[label] = sum(
        (
            coefficient * ct_solution.reduce(power, integral=ct_integral)
            for power, coefficient in ct_parts[label]
        ),
        zero,
    ).together()
    assert (
        ct_reduced[label].replace(d, 4 - 2 * eps).series(eps, 0, -1).to_expression()
        == zero
    ), "CT master coefficient singular at D=4"
    ct_integrated[label] = ct_reduced[label].replace(ct_integral(1), tadpole).replace(
        d, 4 - 2 * eps
    ) * (1 + eps * log_4pi)
    ct_poles[label] = (
        (weight * ct_integrated[label]).series(eps, 0, -1).to_expression().expand()
    )
bare_finite = ct_integrated["bare"].series(eps, 0, 0).to_expression().expand()
ct_sum = sum((ct_poles[label] for label in ct_poles if label != "bare"), zero).expand()
za_one = -flavors * ct_poles["bare"].coefficient(Qpp)
zam_one = (
    -flavors * ct_poles["bare"].coefficient(Qg).replace(s, zero) / (2 * mass_squared)
).together()
# These references are assertions after generation and integration, not inputs.
T = s * Qg - Qpp
assert (ct_poles["bare"] - (4 * mass_squared * Qg - 4 * T / 3) / eps).together() == zero
assert (za_one + 4 * flavors / (3 * eps)).together() == zero
assert (zam_one + 2 * flavors / eps).together() == zero
assert (
    flavors * ct_poles["bare"] - za_one * T + 2 * mass_squared * zam_one * Qg
).together() == zero
expected_ct = xi * mass_squared * Qg * (
    4 / eps**2 + (-4 * log_mass + 4 * log_4pi - 4) / eps
) + 4 * xi * T / (3 * eps)
assert (ct_sum - expected_ct).together() == zero
assert (ct_poles["vertex_0"] - ct_poles["vertex_1"]).expand() == zero
assert (ct_poles["kinetic_0"] - ct_poles["kinetic_1"]).expand() == zero
assert ct_sum.coefficient(eps ** (-3)) == zero
assert ct_sum.replace(xi, zero) == zero

# The auxiliary photon operator in the reference model is i*M*(ZAm²-1)*g.
# Its second-order coefficient includes the square of the derived one-loop term.
loop_sum = (bare_uv + flavors * ct_sum).expand()
za_two = -loop_sum.coefficient(Qpp)
zam_two = -(loop_sum.coefficient(Qg).replace(s, zero) / mass_squared + zam_one**2) / 2
za_two, zam_two = za_two.together(), zam_two.together()
tree_two = -za_two * T + mass_squared * (2 * zam_two + zam_one**2) * Qg
assert (loop_sum + tree_two).together() == zero
assert (za_two + 2 * flavors / eps).together() == zero
assert (
    zam_two - flavors * xi / (2 * eps) + flavors * (2 * flavors + xi) / eps**2
).together() == zero
assert za_two.derivative(xi).expand() == zero
assert za_two.coefficient(eps**-2) == zero
for value in (za_one, zam_one, za_two, zam_two):
    assert value.derivative(log_mass).expand() == zero
    assert value.derivative(log_4pi).expand() == zero
    assert value.derivative(mass_squared).expand() == zero
    assert value.coefficient(eps**-3) == zero
assert (
    loop_sum - flavors * (mass_squared * xi * Qg * (2 / eps**2 - 1 / eps) - 2 * T / eps)
).together() == zero
print(
    "Generated symbolic-gauge two-loop photon UV, four calculated CT insertions, "
    "one-/two-loop photon counterterms and exact tensor cancellation passed",
    {
        "two_loop_targets": len(targets),
        "two_loop_ibp": solution.stats,
        "ct_targets": len(ct_targets),
        "ct_ibp": ct_solution.stats,
    },
)

"""Generated two-loop phi4 self-energy and counterterms with shared IBP reduction.

Reference: https://feyncalc.github.io/FeynCalcExamples/Phi4/TwoLoops/Renormalization-SS
Analytic tadpole/vacuum pole coefficients are reference inputs. A finite
Laporta residual list is not a certification of a minimal master basis.
"""

import json
from math import prod

from symbolica import E, Expression, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

d, k1, k2, mass_squared, p_squared, eps, coupling = S(
    "ibp_phi4::d",
    "ibp_phi4::k1",
    "ibp_phi4::k2",
    "ibp_phi4::M",
    "ibp_phi4::p2",
    "ibp_phi4::eps",
    "ibp_phi4::g",
)
tadpole_squared, vacuum_integral, log_mass, log_4pi = S(
    "ibp_phi4::T2",
    "ibp_phi4::V",
    "ibp_phi4::Lm",
    "ibp_phi4::L4pi",
)
integral = S("ibp_phi4::I")
zero = E("0")
pi = Expression.PI

# Taylor-expand the shifted line, then project its rank-two vacuum numerator.
q, p, scaling, denominator = S(
    "ibp_phi4::q", "ibp_phi4::p", "ibp_phi4::scaling", "ibp_phi4::D3"
)
mink, dot = S("spenso::mink", "spenso::dot")
qc, pc = q(mink(d)), p(mink(d))
shifted_line = 1 / (denominator + 2 * scaling * dot(qc, pc) + scaling**2 * dot(pc, pc))
taylor = shifted_line.series(scaling, 0, 2).to_expression().replace(scaling, E("1"))
taylor = hep.TensorReducer(d).with_integrated_vector(qc).reduce(taylor)
taylor = (
    taylor.replace(dot(qc, qc), denominator + mass_squared)
    .replace(dot(pc, pc), p_squared)
    .expand()
)
expected_taylor = (
    1 / denominator
    + p_squared * (4 / d - 1) / denominator**2
    + 4 * mass_squared * p_squared / (d * denominator**3)
)
assert (taylor - expected_taylor).together() == zero
reference_bare_input = (
    integral(2, 1, 0) / 4
    + taylor.replace_multiple(
        [
            Replacement(denominator**-1, integral(1, 1, 1)),
            Replacement(denominator**-2, integral(1, 1, 2)),
            Replacement(denominator**-3, integral(1, 1, 3)),
        ]
    )
    / 6
)

kinematics = hep.Kinematics(d, momenta=[k1, k2])
family = hep.IntegralFamily(
    [k1, k2],
    [],
    [kinematics.scalar_product(k, k) - mass_squared for k in (k1, k2, k1 + k2)],
    kinematics=kinematics,
)
assert family.is_complete and family.is_independent
# Generate the two bare topologies and carry their native factors through the
# shared graph UV expansion. The supplied integrands are independent checks.
model = hep.Model.phi4()
particle = model.particle("phi")
generated = model.generate_diagrams(
    [particle],
    [particle],
    loops=2,
    max_vertices=2,
    maximum_bridges=0,
    self_energy=None,
    tadpoles=None,
    zero_snails=None,
    numerator_grouping=None,
    progress=None,
)
assert len(generated.diagrams) == 2
K, P, dim = S("gammalooprs::K", "gammalooprs::P", "gammalooprs::dim")
model_mass, model_coupling = S("UFO::mass", "UFO::lam")
index, dimension = S("ibp_phi4::index_", "ibp_phi4::dimension_")
den, edge_, mom_, mass_, quad_ = S(
    "gammalooprs::denom",
    "ibp_phi4::edge_",
    "ibp_phi4::mom_",
    "ibp_phi4::mass_",
    "ibp_phi4::quad_",
)
coordinates = S("ibp_phi4::x1", "ibp_phi4::x2", "ibp_phi4::x3")
external_kinematics = kinematics.with_scalar_product(P(0), P(0), p_squared)
reducer = (
    hep.TensorReducer(d)
    .with_integrated_vector(k1(mink(d)))
    .with_integrated_vector(k2(mink(d)))
)
bare_input = zero
bare_terms = []
raw_checks = []
for diagram in generated.diagrams:
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    factor = (
        diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    raw_family = diagram.integral_family(kinematics=hep.Kinematics(d))
    raw = factor * numerator / (Symbol.I * model_coupling**2)
    for propagator in raw_family.denominators:
        raw /= propagator
    raw = (
        raw.replace(K(0, index), k1(index))
        .replace(K(1, index), k2(index))
        .replace(model_mass**2, mass_squared)
    )
    raw_checks.append(raw)
    expanded = diagram.momentum_basis().route_expression(
        diagram.uv_expansion(model_mass, numerator=numerator).to_expression()
    )
    expanded = (
        expanded.replace(dim, d)
        .replace(K(0, index), k1(index))
        .replace(K(1, index), k2(index))
        .replace(model_mass**2, mass_squared)
    )
    for match in list(expanded.match(den(edge_, mom_, mass_, quad_))):
        values = dict(match)
        expanded = expanded.replace(
            den(values[edge_], values[mom_], values[mass_], values[quad_]),
            family.rewrite_numerator(values[quad_], coordinates),
        )
    scalar = family.rewrite_numerator(
        external_kinematics.apply(reducer.reduce(expanded)), coordinates
    )
    scalar = (scalar * factor / (Symbol.I * model_coupling**2)).together().expand()
    for monomial, coefficient in scalar.coefficient_list(*coordinates):
        powers = [
            -int((monomial.derivative(x) * x / monomial).together())
            for x in coordinates
        ]
        assert monomial == prod(
            x ** (-power) for x, power in zip(coordinates, powers, strict=True)
        )
        assert not coefficient.matches(K(S("ibp_phi4::args___")))
        bare_terms.append((powers, coefficient))
        bare_input += coefficient * integral(*powers)
raw_reference = [
    1 / (4 * family.denominators[0] ** 2 * family.denominators[1]),
    1
    / (
        6
        * family.denominators[0]
        * family.denominators[1]
        * (
            external_kinematics.scalar_product(k1 + k2 + P(0), k1 + k2 + P(0))
            - mass_squared
        )
    ),
]
assert all(
    sum(
        external_kinematics.apply(raw - reference).together() == zero
        for raw in raw_checks
    )
    == 1
    for reference in raw_reference
)
assert (bare_input - reference_bare_input).together() == zero

ibp = hep.IBPFamily(family, name="phi4_vacuum")
identities = ibp.ibp_identities()
assert len(identities) == 4
symmetries = [
    family.mapping_to(family, images)
    for images in (
        [k1, k2],
        [k2, k1],
        [-k1 - k2, k2],
        [k1, -k1 - k2],
        [k2, -k1 - k2],
        [-k1 - k2, k1],
    )
]
assert all(mapping is not None for mapping in symmetries)
assert len({tuple(mapping.denominator_map) for mapping in symmetries}) == 6

targets = sorted({tuple(powers) for powers, coefficient in bare_terms})
assert set(targets) == {(2, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3)}
laporta = ibp.reduce_laporta(targets, max_depth=2)
assert laporta.stats["rows"] > 0
# These are the solver's unresolved integrals at finite search depth. Their
# subsequent identification uses independently verified loop-momentum maps.
canonical = {
    tuple(powers): min(tuple(mapping.map_powers(powers)) for mapping in symmetries)
    for powers in laporta.residuals
}
reference_basis = {(0, 1, 1): tadpole_squared, (1, 1, 1): vacuum_integral}
raw_reductions = {tuple(target): laporta.reduce(target) for target in targets}
reduced = {}
for target, terms in raw_reductions.items():
    assert all(tuple(powers) in canonical for powers, coefficient in terms)
    assert all(
        canonical[tuple(powers)] in reference_basis for powers, coefficient in terms
    )
    reduced[target] = sum(
        (
            coefficient * reference_basis[canonical[tuple(powers)]]
            for powers, coefficient in terms
        ),
        zero,
    ).together()
assert (
    reduced[2, 1, 0] - (d - 2) * tadpole_squared / (2 * mass_squared)
).together() == zero
assert (
    reduced[1, 1, 2] - (d - 3) * vacuum_integral / (3 * mass_squared)
).together() == zero
assert (
    reduced[1, 1, 3]
    - (d - 8) * (d - 3) * vacuum_integral / (18 * mass_squared**2)
    - (d - 2) ** 2 * tadpole_squared / (12 * mass_squared**3)
).together() == zero
bare_reduced = bare_input.replace_multiple(
    [Replacement(integral(*target), value) for target, value in reduced.items()]
).together()

# Analytic master coefficients are reference inputs, not computed by IBP.
tadpole_poles = mass_squared * (1 / eps + 1 - log_mass)
tadpole_squared_poles = mass_squared**2 * (1 / eps**2 + 2 * (1 - log_mass) / eps)
vacuum_poles = mass_squared * (E("3/2") / eps**2 + (E("9/2") - 3 * log_mass) / eps)
bare_uv = (
    (
        -(1 + 2 * eps * log_4pi)
        * bare_reduced.replace(d, 4 - 2 * eps)
        .replace(tadpole_squared, tadpole_squared_poles)
        .replace(vacuum_integral, vacuum_poles)
    )
    .series(eps, 0, -1)
    .to_expression()
    .expand()
)
expected_bare = (
    -mass_squared / (2 * eps**2)
    + (mass_squared * (log_mass - 1 - log_4pi) + p_squared / 24) / eps
)
assert (bare_uv - expected_bare).together() == zero

tadpole_kin = hep.Kinematics(d, momenta=[k1])
tadpole_family = hep.IntegralFamily(
    [k1],
    [],
    [tadpole_kin.scalar_product(k1, k1) - mass_squared],
    kinematics=tadpole_kin,
)
parametric = hep.IBPFamily(tadpole_family, name="phi4_tadpole").solve_parametric(
    [True], max_depth=1
)
assert parametric.rules
tadpole_terms = parametric.reduce([2])
assert len(tadpole_terms) == 1 and tadpole_terms[0][0] == [1]
tadpole_coefficient = tadpole_terms[0][1]
assert (tadpole_coefficient - (d - 2) / (2 * mass_squared)).together() == zero
zg1, zm1 = 3 / (2 * eps), 1 / (2 * eps)
# Build local counterterm vertices from the bare Lagrangian factors. Keep the
# one-loop field constant symbolic through graph generation, although its
# reference value vanishes in phi4. h counts powers of g/(16*pi^2).
h, field1, mass1, vertex1, field2, mass2 = S(
    "ibp_phi4::h",
    "ibp_phi4::field1",
    "ibp_phi4::mass1",
    "ibp_phi4::vertex1",
    "ibp_phi4::field2",
    "ibp_phi4::mass2",
)
Zfield = 1 + h * field1 + h**2 * field2
Zmass = 1 + h * mass1 + h**2 * mass2
Zvertex = 1 + h * vertex1
specification = json.loads(model.to_json())
specification["orders"].append({"name": "CT", "expansion_order": 2, "hierarchy": 1})
for label, valence, lorentz, bare_factor in [
    ("kinetic", 2, "P(dummy(1),1)*P(dummy(1),1)", Symbol.I * (Zfield - 1)),
    ("mass", 2, "1", -Symbol.I * model_mass**2 * (Zfield * Zmass - 1)),
    ("quartic", 4, "1", -Symbol.I * model_coupling * (Zvertex * Zfield**2 - 1)),
]:
    specification["lorentz_structures"].append(
        {"name": "CT_L_" + label, "spins": [1] * valence, "structure": lorentz}
    )
    for order in [1, 2] if valence == 2 else [1]:
        name = f"CT_{label}_{order}"
        specification["couplings"].append(
            {
                "name": name,
                "expression": repr(
                    bare_factor.expand().coefficient(h**order) * h**order
                ),
                "orders": [["SCALAR", 1 if valence == 4 else 0], ["CT", order]],
                "value": None,
            }
        )
        specification["vertex_rules"].append(
            {
                "name": name,
                "particles": [particle.name] * valence,
                "color_structures": ["1"],
                "lorentz_structures": ["CT_L_" + label],
                "couplings": [[name]],
            }
        )
ct_model = hep.Model.from_json(json.dumps(specification))
ct_diagrams = {}
ct_integral, ct_coordinate = S("ibp_phi4::J", "ibp_phi4::y")
ct_input, tree_ct = zero, zero
for loops, ct_order in [(1, 1), (0, 2)]:
    result = ct_model.generate_diagrams(
        [particle.name],
        [particle.name],
        loops=loops,
        max_vertices=2 if loops else 1,
        coupling_orders={"CT": ct_order},
        maximum_bridges=0,
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(result.diagrams) == (3 if loops else 2)
    ct_diagrams[loops] = result.diagrams
    for diagram in result.diagrams:
        numerator = ct_model.expand_couplings(
            diagram.numerator_expression().to_expression()
        )
        numerator = TensorExpression(numerator.expand()).to_dots().to_expression()
        if loops:
            numerator = diagram.uv_expansion(
                model_mass, numerator=numerator
            ).to_expression()
        numerator = diagram.momentum_basis().route_expression(numerator)
        numerator = (
            numerator.replace(mink(dimension), mink(d))
            .replace(K(0, index), k1(index))
            .replace(model_mass**2, mass_squared)
        )
        for match in list(numerator.match(den(edge_, mom_, mass_, quad_))):
            values = dict(match)
            numerator = numerator.replace(
                den(values[edge_], values[mom_], values[mass_], values[quad_]),
                tadpole_family.rewrite_numerator(values[quad_], [ct_coordinate]),
            )
        factor = (
            diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
        )
        numerator = external_kinematics.apply(numerator) * factor
        if loops:
            scalar = (
                (
                    tadpole_family.rewrite_numerator(numerator, [ct_coordinate])
                    / (h * model_coupling)
                )
                .together()
                .expand()
            )
            for monomial, coefficient in scalar.coefficient_list(ct_coordinate):
                power = -int(
                    (
                        monomial.derivative(ct_coordinate) * ct_coordinate / monomial
                    ).together()
                )
                assert monomial == ct_coordinate ** (-power)
                ct_input += coefficient * ct_integral(power)
        else:
            tree_ct += numerator / (Symbol.I * h**2)
assert (
    ct_input
    - (vertex1 + field1) * ct_integral(1) / 2
    - mass_squared * mass1 * ct_integral(2) / 2
).together() == zero
assert (
    tree_ct - p_squared * field2 + mass_squared * (mass2 + field2 + field1 * mass1)
).together() == zero
ct_reduced = ct_input.replace(ct_integral(2), tadpole_coefficient * ct_integral(1))
ct_reduced = ct_reduced.replace(field1, zero).replace(vertex1, zg1).replace(mass1, zm1)

counterterm_uv = (
    (
        (1 + eps * log_4pi)
        * ct_reduced.replace(d, 4 - 2 * eps).replace(ct_integral(1), tadpole_poles)
    )
    .series(eps, 0, -1)
    .to_expression()
    .expand()
)
expected_counterterm = (
    mass_squared / eps**2 + mass_squared * (E("3/4") - log_mass + log_4pi) / eps
)
assert (counterterm_uv - expected_counterterm).together() == zero
loop_uv = (bare_uv + counterterm_uv).expand()
assert (
    loop_uv
    - mass_squared / (2 * eps**2)
    + mass_squared / (4 * eps)
    - p_squared / (24 * eps)
).together() == zero
tree_ct = tree_ct.replace(field1, zero).expand()
ct_rows = [
    tree_ct.coefficient(p_squared),
    (tree_ct.replace(p_squared, zero) / mass_squared).together().expand(),
]
ct_matrix = Matrix.from_linear(
    2, 2, [row.coefficient(unknown) for row in ct_rows for unknown in (field2, mass2)]
)
rhs = Matrix.vec(
    [-loop_uv.coefficient(p_squared), -loop_uv.replace(p_squared, zero) / mass_squared]
)
solution = ct_matrix.solve(rhs)
zphi2, zm2 = (solution[row, 0].to_expression().expand() for row in range(2))
assert (loop_uv + tree_ct.replace(field2, zphi2).replace(mass2, zm2)).together() == zero
assert (loop_uv + p_squared * zphi2 - mass_squared * (zm2 + zphi2)).together() == zero
loop_coupling = coupling / (16 * pi**2)
zphi = (1 + loop_coupling**2 * zphi2).expand()
zm = (1 + loop_coupling * zm1 + loop_coupling**2 * zm2).expand()
assert (zphi - 1 + coupling**2 / (6144 * pi**4 * eps)).together() == zero
assert (
    zm
    - 1
    - coupling / (32 * pi**2 * eps)
    - coupling**2 / (512 * pi**4 * eps**2)
    + 5 * coupling**2 / (6144 * pi**4 * eps)
).together() == zero
print(
    "Generated phi4 diagrams and counterterms, native UV/tensor expansion, IBP and renormalization constants passed",
    laporta.stats,
)

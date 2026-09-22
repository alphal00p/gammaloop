"""Shared IBP reduction of the FeynCalc two-loop phi4 self-energy.

Reference: https://feyncalc.github.io/FeynCalcExamples/Phi4/TwoLoops/Renormalization-SS
Analytic tadpole/vacuum pole coefficients are reference inputs. A finite
Laporta residual list is not a certification of a minimal master basis.
"""

from symbolica import E, Expression, Replacement, S
from symbolica.community import hep

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
bare_input = (
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

targets = [[2, 1, 0], [1, 1, 1], [1, 1, 2], [1, 1, 3]]
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
counterterm_uv = (
    (
        (1 + eps * log_4pi)
        * (zg1 + mass_squared * zm1 * tadpole_coefficient.replace(d, 4 - 2 * eps))
        * tadpole_poles
        / 2
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
zphi2 = -loop_uv.coefficient(p_squared).expand()
zm2 = (loop_uv.replace(p_squared, zero) / mass_squared - zphi2).expand()
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
    "Phi4 Laporta and parametric IBP, Taylor projection, UV poles and renormalization constants passed",
    laporta.stats,
)

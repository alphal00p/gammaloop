"""Unequal-mass bubble IBP checked by identities, UV poles and quadrature.

D1=k^2-a and D2=(k-p)^2-b use squared masses a,b and p^2=s.
The finite Laporta residual basis is evaluated with the shared OneLOop A0/B0
convention. Independent Feynman-parameter quadrature checks all three raised
integrals at spacelike and timelike points below the two-particle threshold.
"""

import numpy as np
from symbolica import E, Replacement, S
from symbolica.community import hep
from symbolica.community.hep import oneloop

d, k, p, s, mass_a, mass_b, eps = S(
    "ibp_bubble::d",
    "ibp_bubble::k",
    "ibp_bubble::p",
    "ibp_bubble::s",
    "ibp_bubble::a",
    "ibp_bubble::b",
    "ibp_bubble::eps",
)
tadpole_a, tadpole_b, bubble = S(
    "ibp_bubble::T_a", "ibp_bubble::T_b", "ibp_bubble::B11"
)
finite_a, finite_b, finite_bubble = S(
    "ibp_bubble::A_finite", "ibp_bubble::B_finite", "ibp_bubble::C_finite"
)
kinematics = hep.Kinematics(d, momenta=[k, p]).with_scalar_product(p, p, s)
family = hep.IntegralFamily(
    [k],
    [p],
    [
        kinematics.scalar_product(k, k) - mass_a,
        kinematics.scalar_product(k - p, k - p) - mass_b,
    ],
    kinematics=kinematics,
)
assert family.is_complete and family.is_independent
targets = [[2, 1], [1, 2], [2, 2]]
solution = hep.IBPFamily(family, name="unequal_mass_bubble").reduce_laporta(
    targets, max_depth=2
)
assert solution.stats["rows"] > 0
assert {tuple(powers) for powers in solution.residuals} == {(1, 0), (0, 1), (1, 1)}
reference_basis = {(1, 0): tadpole_a, (0, 1): tadpole_b, (1, 1): bubble}
integral = S("ibp_bubble::I")
basis_replacements = [
    Replacement(integral(*powers), value) for powers, value in reference_basis.items()
]

# Independently solve the two total-derivative identities for I(2,1).
# The discriminant is nonzero at the numerical points below.
discriminant = (
    s**2 + mass_a**2 + mass_b**2 - 2 * s * mass_a - 2 * s * mass_b - 2 * mass_a * mass_b
)
expected_21 = (
    (d - 2) / discriminant * tadpole_b
    + (d - 2) * (s - mass_a - mass_b) / (2 * mass_a * discriminant) * tadpole_a
    - (d - 3) * (s - mass_a + mass_b) / discriminant * bubble
)
exchange = [
    Replacement(mass_a, mass_b),
    Replacement(mass_b, mass_a),
    Replacement(tadpole_a, tadpole_b),
    Replacement(tadpole_b, tadpole_a),
]
expected_12 = expected_21.replace_multiple(exchange)
# Differentiating I(2,1) with respect to b raises its second line. For the
# residuals, dT_b/db=(d-2)T_b/(2b) and dB11/db=I(1,2).
expected_22 = (
    expected_21.derivative(mass_b)
    + expected_21.coefficient(tadpole_b) * (d - 2) * tadpole_b / (2 * mass_b)
    + expected_21.coefficient(bubble) * expected_12
)
expected = {(2, 1): expected_21, (1, 2): expected_12, (2, 2): expected_22}

# OneLOop returns Laurent coefficients [finite, 1/eps, 1/eps^2].
# Recombine poles before d=4-2eps: O(eps) reduction coefficients multiply
# divergent residual integrals and contribute to the finite result.
poles = [
    Replacement(tadpole_a, finite_a + mass_a / eps),
    Replacement(tadpole_b, finite_b + mass_b / eps),
    Replacement(bubble, finite_bubble + 1 / eps),
]
nodes, weights = np.polynomial.legendre.leggauss(96)
nodes, weights = (nodes + 1) / 2, weights / 2
points = [(-1.0, 2.0, 3.0), (2.0, 2.0, 3.0), (5.0, 1.0, 4.0), (8.0, 1.0, 4.0)]

for target in targets:
    terms = solution.reduce(target)
    assert all(tuple(powers) in reference_basis for powers, coefficient in terms)
    expression = solution.reduce(target, integral=integral)
    tuple_expression = sum(
        (coefficient * integral(*powers) for powers, coefficient in terms), E("0")
    )
    assert (expression - tuple_expression).together() == E("0")
    reduced = expression.replace_multiple(basis_replacements)
    assert (reduced - expected[tuple(target)]).together() == E("0"), target
    if target == [2, 2]:
        assert (reduced - reduced.replace_multiple(exchange)).together() == E("0")
    series = (
        reduced.replace_multiple(poles)
        .replace(d, 4 - 2 * eps)
        .series(eps, 0, 0)
        .to_expression()
        .expand()
    )
    assert series.coefficient(eps**-1).together() == E("0")
    assert series.coefficient(eps**-2).together() == E("0")
    finite = dict(series.coefficient_list(eps))[E("1")].together()
    for invariant, first_mass, second_mass in points:
        a0 = oneloop.A0(first_mass, 1.0)
        a1 = oneloop.A0(second_mass, 1.0)
        b0 = oneloop.B0(invariant, first_mass, second_mass, 1.0)
        assert abs(complex(a0[1]) - first_mass) < 1e-12
        assert abs(complex(a1[1]) - second_mass) < 1e-12
        assert abs(complex(b0[1]) - 1) < 1e-12
        numeric = complex(
            finite.evaluate(
                {
                    s: invariant,
                    mass_a: first_mass,
                    mass_b: second_mass,
                    finite_a: complex(a0[0]),
                    finite_b: complex(a1[0]),
                    finite_bubble: complex(b0[0]),
                }
            )
        )
        denominator = (
            nodes * first_mass
            + (1 - nodes) * second_mass
            - invariant * nodes * (1 - nodes)
        )
        assert np.all(denominator > 0)
        # The normalized Minkowski parameter formula has sign (-1)^(n1+n2):
        # negative for I(2,1), I(1,2), and positive for I(2,2).
        if target == [2, 1]:
            integrand = -nodes / denominator
        elif target == [1, 2]:
            integrand = -(1 - nodes) / denominator
        else:
            integrand = nodes * (1 - nodes) / denominator**2
        independent = float(np.dot(weights, integrand))
        assert abs(numeric - independent) < 2e-11, (
            target,
            invariant,
            numeric,
            independent,
        )
    print(
        f"Bubble {target}: analytic IBP, exact UV cancellation and four quadrature checks passed",
        flush=True,
    )

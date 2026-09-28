"""Exact propagator decompositions through the installed family API.

Reference algebra: https://feyncalc.github.io/FeynCalcBookDev/ApartFF.html
Momentum shifts and scaleless-integral removal are deliberately separate steps.
"""

import importlib
import sys

from symbolica import E, S

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)

k, p, s = S("apart_check::k", "apart_check::p", "s")
kin = fk.Kinematics(momenta=[k, p]).with_scalar_product(p, p, s)
kk = kin.scalar_product(k, k)
kp = kin.scalar_product(k, p)

for denominators, powers in [
    (
        [kk, kin.scalar_product(k - p, k - p), kin.scalar_product(k + p, k + p)],
        [1, 1, 1],
    ),
    ([kk, kp, kk + kp], [2, 1, 2]),
    ([kk, kk - s, kp, kk + kp], [2, 3, -1, 1]),
]:
    family = fk.IntegralFamily([k], [p], denominators, kinematics=kin)
    terms = family.partial_fraction(powers)
    expected = E("1")
    for denominator, power in zip(denominators, powers):
        expected *= denominator ** (-power)
    reconstructed = E("0")
    for coefficient, term_powers in terms:
        term = coefficient
        active = []
        for denominator, power in zip(denominators, term_powers):
            term *= denominator ** (-power)
            if power > 0:
                active.append(denominator)
        assert fk.IntegralFamily([k], [p], active, kinematics=kin).is_independent
        reconstructed += term
    assert (reconstructed - expected).together() == E("0")

print("Partial fractions preserve rational functions and independent supports")

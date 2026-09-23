"""Loop-family completion and exact numerator mappings in the installed host."""

from symbolica import E, S
from symbolica.community import hep as fk

D, k, q, p, s, m1, m2, d1, d2 = S(
    "family::D", "family::k", "family::q", "family::p", "s", "m1sq", "m2sq", "d1", "d2"
)
kin = fk.Kinematics(D, momenta=[k, q, p]).with_scalar_product(p, p, s)
kk = kin.scalar_product(k, k)
kp = kin.scalar_product(k, p)
bubble = fk.IntegralFamily(
    [k], [p], [kk - m1, kin.scalar_product(k + p, k + p) - m2], kinematics=kin
)
assert bubble.rank == 2 and bubble.is_complete and bubble.is_independent
reduced = bubble.rewrite_numerator(kp**2, [d1, d2])
assert (reduced - (d2 - d1 - s + m2 - m1) ** 2 / 4).expand() == E("0")

family = fk.IntegralFamily(
    [k, q],
    [p],
    [kin.scalar_product(v, v) for v in [k, q, k - p, q - p]],
    kinematics=kin,
)
assert family.rank == 4 and not family.is_complete
completed = family.complete()
assert completed.rank == 5 and completed.is_complete
assert completed.denominators[:4] == family.denominators
assert completed.denominators[4] == kin.scalar_product(k, q)
preferred = kin.scalar_product(k - q, k - q)
selected = family.complete(candidates=[family.denominators[0], preferred])
assert selected.denominators == [*family.denominators, preferred]
assert selected.complete(candidates=[preferred]).denominators == selected.denominators
assert (
    family.complete(candidates=[family.denominators[0]]).denominators
    == completed.denominators
)
for candidate_family in (family, selected):
    try:
        candidate_family.complete(candidates=[preferred, kk**2])
    except fk.IntegralFamilyError:
        pass
    else:
        raise AssertionError("nonlinear completion candidate was silently accepted")

dependent = fk.IntegralFamily([k], [], [kk, kk - m1], kinematics=kin)
assert dependent.is_complete and not dependent.is_independent
try:
    dependent.complete()
except fk.IntegralFamilyError:
    pass
else:
    raise AssertionError("dependent propagators were accepted for completion")
print("Integral-family rank, completion and numerator mapping checks passed")

# Cofactors define U/F even where completing the square has no inverse.
x, y, delta = S("family::x", "family::y", "family::delta")
singular = fk.IntegralFamily(
    [k, q],
    [p],
    [kin.scalar_product(k + q, k + q) - m1, kin.scalar_product(k - q, p) + delta],
    kinematics=kin,
)
U, F = singular.symanzik([x, y])
assert U == E("0")
assert (F - s * x * y**2).expand() == E("0")
for operation in (
    lambda: singular.scaleless_scaling([x, y]),
    lambda: singular.parametric_mapping(singular, [x, y]),
):
    try:
        operation()
    except fk.IntegralFamilyError:
        pass
    else:
        raise AssertionError(
            "singular U/F polynomials were accepted as a scaling or mapping proof"
        )
print("Singular Symanzik polynomials and proof-domain checks passed")

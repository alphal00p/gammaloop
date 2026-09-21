"""Loop-family completion and exact numerator mappings in the installed host."""

from symbolica import E, S
from symbolica.community import feynkit as fk

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

dependent = fk.IntegralFamily([k], [], [kk, kk - m1], kinematics=kin)
assert dependent.is_complete and not dependent.is_independent
try:
    dependent.complete()
except fk.IntegralFamilyError:
    pass
else:
    raise AssertionError("dependent propagators were accepted for completion")
print("Integral-family rank, completion and numerator mapping checks passed")

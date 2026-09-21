"""Momentum combinations and Mandelstam conservation in the installed host."""

from symbolica import E, S
from symbolica.community import feynkit as fk

p1, p2, p3, p4, s, t, a = S("p1", "p2", "p3", "p4", "s", "t", "a")
masses = list(S("m1sq", "m2sq", "m3sq", "m4sq"))
u = sum(masses) - s - t
kin = fk.Kinematics.mandelstam([p1, p2, p3, p4], masses, [s, t, u])
assert kin.scalar_product(p1 + p2, p1 + p2) == s
assert kin.scalar_product(p1 - p3, p1 - p3) == t
for p in [p1, p2, p3, p4]:
    assert kin.scalar_product(p, p1 + p2 - p3 - p4) == E("0")
assert (kin.scalar_product(a * (p1 + p2), p3 + p4) - a * s).expand() == E("0")

free = fk.Kinematics(momenta=[p1, p2])
assert free.scalar_product(p1 + p2, p1 - p2) == (
    free.scalar_product(p1, p1) - free.scalar_product(p2, p2)
)
for invalid in [p1 * p2, p1**2, p1 / p2, p1 + 1]:
    try:
        free.scalar_product(invalid, p1)
    except fk.KinematicsError:
        pass
    else:
        raise AssertionError(f"accepted nonlinear or nonvector input: {invalid}")
print("Symbolic momentum and Mandelstam conservation checks passed")

"""Momentum combinations and Mandelstam conservation in the installed host."""

from symbolica import E, Expression, S
from symbolica.community import hep as fk
from symbolica.community.spenso import Representation, TensorExpression, TensorName

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

# Kinematic assumptions can annihilate an open tensor. Its ordered ports must
# survive even though the zero atom alone has no recoverable tensor structure.
lorentz = Representation.mink(4)
mu, nu = lorentz("kin_mu"), lorentz("kin_nu")
vector = TensorName.vector("kin_vector")
ordered = vector(nu) * vector(mu)
dot = free.scalar_product(p1, p2)
assumed = free.with_scalar_product(p1, p2, s)
for tensor in (
    dot * ordered,
    (dot - s) * ordered,
    0 * ordered,
    TensorExpression(dot),
):
    result = assumed.apply(tensor)
    assert isinstance(result, TensorExpression)
    assert result.structure.slots == tensor.structure.slots
    plain = assumed.apply(tensor.to_expression())
    assert isinstance(plain, Expression) and not isinstance(plain, TensorExpression)
    assert result.to_expression() == plain
    assert assumed.apply(result) == result
assert assumed.apply((dot - s) * ordered).rank == 2
assert assumed.apply((dot - s) * ordered) == 0
assert assumed.apply(dot * ordered) == s * ordered
for scalar in (E("3"), 3, 2.5):
    result = assumed.apply(scalar)
    assert isinstance(result, Expression) and not isinstance(result, TensorExpression)
    assert result == scalar
print("Symbolic momentum, Mandelstam conservation and typed kinematics checks passed")

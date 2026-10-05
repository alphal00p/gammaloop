"""Shared numerical external states exposed through the installed HEPKit host."""

import math
from symbolica.community import hepkit as hep

p = hep.FourMomentum(2.0, 0.0, 0.0, -2.0)
q = 1 / math.sqrt(2)
epsilon = p.wavefunction("epsilon", hep.Helicity.PLUS)
assert epsilon.kind == "epsilon"
assert len(epsilon) == 4
for actual, expected in zip(epsilon.components, [0j, -q, 1j * q, 0j]):
    assert abs(actual - expected) < 1e-14
assert epsilon.bar().kind == "epsilon_bar"
assert epsilon.bar().components == [x.conjugate() for x in epsilon.components]
assert epsilon.bar().bar() == epsilon

for kind, helicity, expected in [
    ("u", hep.Helicity.PLUS, [0, 0, 0, -2]),
    ("u", hep.Helicity.MINUS, [2, 0, 0, 0]),
    ("v", hep.Helicity.PLUS, [-2, 0, 0, 0]),
    ("v", hep.Helicity.MINUS, [0, 0, 0, 2]),
]:
    state = p.wavefunction(kind, helicity)
    assert state.components == expected
    adjoint = state.bar()
    assert adjoint.kind == kind + "_bar"
    assert adjoint.components == [complex(x).conjugate() for x in expected[2:] + expected[:2]]
    assert p.wavefunction(kind + "_bar", helicity) == adjoint
    assert adjoint.bar() == state

massive = hep.FourMomentum(5.0, 0.0, 0.0, 3.0)
assert massive.wavefunction("epsilon", hep.Helicity.ZERO).components == [0.75, 0, 0, 1.25]
scalar = massive.wavefunction("scalar", hep.Helicity.ZERO)
assert scalar.kind == "scalar" and scalar.components == [1] and len(scalar) == 1
assert scalar.bar() == scalar

for momentum, kind, helicity in [
    (p, "epsilon", hep.Helicity.ZERO),
    (p, "u", hep.Helicity.ZERO),
    (p, "scalar", hep.Helicity.PLUS),
    (p, "unknown", hep.Helicity.PLUS),
    (hep.FourMomentum(float("nan"), 0, 0, 1), "epsilon", hep.Helicity.PLUS),
]:
    try:
        momentum.wavefunction(kind, helicity)
    except hep.KinematicsError:
        pass
    else:
        raise AssertionError(f"accepted invalid {kind} state")

# Returned Python containers cannot modify the native immutable state.
copy = epsilon.components
copy[0] = 123j
assert epsilon.components[0] == 0j
print("Shared scalar, vector, spinor and adjoint wavefunctions passed")

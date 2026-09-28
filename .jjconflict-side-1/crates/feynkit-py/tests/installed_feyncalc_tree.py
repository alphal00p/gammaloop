"""QED spin sums through the installed FeynKit, Idenso and Symbolica APIs.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-MuAmu
This tests the squared current contraction, not diagram generation or phase space.
"""

import importlib
import sys
from pathlib import Path

from symbolica import E, S
from symbolica.community.spenso import (
    Representation,
    TensorExpression,
    TensorName,
    chain,
)

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
p1, p2, k1, k2 = (
    TensorName.vector(name).to_expression() for name in ("p1", "p2", "k1", "k2")
)
s, t, u, e = S("s", "t", "u", "e")
me, mm = S("UFO::Me", "UFO::MM")
gamma, bis, mink = S("spenso::gamma", "spenso::bis", "spenso::mink")
mu, nu = S("mu", "nu")
i0, i1, i2, i3, j0, j1, j2, j3 = S("i0", "i1", "i2", "i3", "j0", "j1", "j2", "j3")

incoming = (
    model.particle_by_pdg(-11).spin_sum(p2, i0, i1, average=True)
    * gamma(bis(4, i1), bis(4, i2), mink(4, mu))
    * model.particle_by_pdg(11).spin_sum(p1, i2, i3, average=True)
    * gamma(bis(4, i3), bis(4, i0), mink(4, nu))
)
outgoing = (
    model.particle_by_pdg(13).spin_sum(k1, j0, j1)
    * gamma(bis(4, j1), bis(4, j2), mink(4, mu))
    * model.particle_by_pdg(-13).spin_sum(k2, j2, j3)
    * gamma(bis(4, j3), bis(4, j0), mink(4, nu))
)
kin = fk.Kinematics.mandelstam(
    [p1, p2, k1, k2], [me**2, me**2, mm**2, mm**2], [s, t, u]
)
contracted = (
    (
        TensorExpression(incoming).simplify_gamma().to_expression()
        * TensorExpression(outgoing).simplify_gamma().to_expression()
    )
    .contract()
    .to_dots()
    .to_expression()
    .to_expression()
)
squared = kin.apply(contracted) * e**4 / s**2
expected = (
    2
    * e**4
    / s**2
    * (
        2 * me**2 * (2 * mm**2 + s - t - u)
        + 2 * me**4
        + 2 * mm**4
        + 2 * mm**2 * (s - t - u)
        + t**2
        + u**2
    )
)
assert (squared - expected).together() == E("0"), squared
massless = squared.replace(me, E("0")).replace(mm, E("0"))
assert (massless - 2 * e**4 * (t**2 + u**2) / s**2).together() == E("0")

# Independent contexts retain their own assumptions.
alternative = fk.Kinematics().with_scalar_product(p1, p2, E("7"))
assert alternative.scalar_product(p1, p2) == E("7")
assert kin.scalar_product(p1, p2) == s / 2 - me**2
print("FeynCalc QED massive and massless squared-current benchmarks passed")


# Right-chiral currents: P_L gamma^mu P_R = gamma^mu P_R. For massless
# external states these select right-handed particles and left-handed antiparticles.
# The typed gamma pattern leaves the momentum slashes in the spin sums intact.
a, b, ell = S("polarized::a_", "polarized::b_", "polarized::ell_")
pattern = TensorExpression.gamma(4)(a, b, ell).to_expression()
right_current = chain(
    Representation.bis(4)(a),
    Representation.bis(4)(b),
    TensorExpression.gamma(4)(a, "middle", ell),
    TensorExpression.projp(4)("middle", b),
).to_expression()
projected = (incoming * outgoing).replace(pattern, right_current)
traces = TensorExpression(projected).simplify_gamma()
contracted = (
    traces.simplify_epsilon().contract().to_dots().to_expression().to_expression()
)
# Undo the two initial spin averages: each projected state is specified.
polarized = 4 * kin.apply(contracted) * e**4 / s**2
assert (polarized - 4 * e**4 * (me**2 + mm**2 - u) ** 2 / s**2).together() == E("0")
polarized_massless = polarized.replace(me, E("0")).replace(mm, E("0"))
assert (polarized_massless - 4 * e**4 * u**2 / s**2).together() == E("0")
print("FeynCalc QED polarized massive and massless squared-current benchmarks passed")

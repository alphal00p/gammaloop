"""Independent review: Grassmann signs, derivatives, conjugation and index spaces, on small cases worked out by hand.

Each expected result is derived in the docstring from the conventions of the report (02_method sec:m-objects,
03b sec:st-signs: left derivative ``F[X + dX] = F + int dX dL F/dX``, right derivative ``F + int (dR F/dX) dX``,
graded left derivations, vertex function = coefficient of ``theta_1 ... theta_n``, hermitian conjugation reverses
the order of the factors, ``Psibar^dagger -> Psi`` with the Dirac adjoint, ``omega^dagger = omega``,
``omegabar^dagger = -omegabar``, ``zeta^dagger = -zeta``).
"""

from __future__ import annotations

import pytest
from symbolica import Expression

from hepkit_renorm import builtin_model
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST
from hepkit_renorm.functional import Density, LocalFunctional, Occurrence, standard_legs


@pytest.fixture(scope="module")
def sm():
    return builtin_model("sm_one_family")


@pytest.fixture(scope="module")
def simp(sm):
    return sm.simplifier()


def labels(sm, name):
    return sm.slot_labels(name)


def mono(sm, coeff, *occs):
    return Density.monomial(sm.fields, coeff, *occs)


def test_reordering_and_antisymmetry(sm):
    """omega^a omega^b = - omega^b omega^a; delta^{ab} omega^a omega^b = 0; Qbar Q = - Q Qbar; Qbar H = H Qbar."""
    f = sm.fields
    a, b = labels(sm, "gh_G")[0], labels(sm, "gh_G")[0]
    ab = mono(sm, 1, Occurrence("gh_G", (a,)), Occurrence("gh_G", (b,)))
    ba = mono(sm, 1, Occurrence("gh_G", (b,)), Occurrence("gh_G", (a,)))
    assert ab == -ba and not ab.is_zero()
    assert mono(sm, tn.delta(a, b), Occurrence("gh_G", (a,)), Occurrence("gh_G", (b,))).is_zero()
    c = labels(sm, "gh_G")[0]
    fabc = mono(sm, tn.f(a, b, c), Occurrence("gh_G", (b,)), Occurrence("gh_G", (c,)))
    assert not fabc.is_zero()
    q = labels(sm, "Q")
    qb = labels(sm, "Qbar")
    one = tn.delta(qb[0], q[0]) * tn.delta(qb[1], q[1]) * tn.delta(qb[2], q[2])
    assert mono(sm, one, Occurrence("Qbar", qb), Occurrence("Q", q)) == \
        -mono(sm, one, Occurrence("Q", q), Occurrence("Qbar", qb))
    h = labels(sm, "H")
    assert mono(sm, 1, Occurrence("Qbar", qb), Occurrence("H", h)) == mono(sm, 1, Occurrence("H", h),
                                                                           Occurrence("Qbar", qb))
    del f


def test_graded_product(sm):
    """For odd densities X, Y with disjoint labels X * Y = - Y * X; for an even one the factors commute."""
    x = mono(sm, 1, Occurrence("gh_B"))
    y = mono(sm, 1, Occurrence("Ebar", labels(sm, "Ebar")))
    h = mono(sm, 1, Occurrence("H", labels(sm, "H")))
    assert x * y == -(y * x) and not (x * y).is_zero()
    assert x * h == h * x


def test_left_and_right_derivatives(sm, simp):
    """F = int Qbar_s M_st Q_t, M = gammabar^mu P_L G_mu-free part replaced by a constant matrix (here PL):
    dL F/dQbar = M Q; dR F/dQbar = - M Q; dL F/dQ = - Qbar M; dR F/dQ = + Qbar M.
    With a derivative, F = int Qbar gamma^mu d_mu Q: dL F/dQ_t = (d_mu Qbar_s) gamma^mu_st (integration by parts:
    Qbar dQ = -(dQ) Qbar ~ + Q d(Qbar))."""
    qb, q = labels(sm, "Qbar"), labels(sm, "Q")
    deltas = tn.delta(qb[1], q[1]) * tn.delta(qb[2], q[2])
    F = LocalFunctional.integrate(mono(sm, tn.PL(qb[0], q[0]) * deltas, Occurrence("Qbar", qb), Occurrence("Q", q)))
    ext = labels(sm, "Qbar")
    mq = mono(sm, tn.PL(ext[0], q[0]) * tn.delta(ext[1], q[1]) * tn.delta(ext[2], q[2]), Occurrence("Q", q))
    assert F.functional_derivative("Qbar", ext, "left").simplify(simp) == mq.simplify(simp)
    assert F.functional_derivative("Qbar", ext, "right").simplify(simp) == (-mq).simplify(simp)
    ext = labels(sm, "Q")
    qbm = mono(sm, tn.PL(qb[0], ext[0]) * tn.delta(qb[1], ext[1]) * tn.delta(qb[2], ext[2]), Occurrence("Qbar", qb))
    assert F.functional_derivative("Q", ext, "left").simplify(simp) == (-qbm).simplify(simp)
    assert F.functional_derivative("Q", ext, "right").simplify(simp) == qbm.simplify(simp)
    mu = tn.fresh(tn.LORENTZ)
    K = LocalFunctional.integrate(mono(sm, tn.gamma(mu, qb[0], q[0]) * deltas, Occurrence("Qbar", qb),
                                       Occurrence("Q", q, (mu,))))
    ext = labels(sm, "Q")
    nu = tn.fresh(tn.LORENTZ)
    exp = mono(sm, tn.gamma(nu, qb[0], ext[0]) * tn.delta(qb[1], ext[1]) * tn.delta(qb[2], ext[2]),
               Occurrence("Qbar", qb, (nu,)))
    assert K.functional_derivative("Q", ext, "left").simplify(simp) == exp.simplify(simp)


def test_odd_source_right_derivative(sm, simp):
    """L_ext contains (s H^dag) sigma with s H^dag and sigma odd: dL/dsigma = - s H^dag, dR/dsigma = + s H^dag
    (03b: the right derivative in the H^dag term of eq:ST)."""
    br = BRST(sm)
    lab = labels(sm, "src_Hdag")
    shd = br.variation("Hdag", lab).simplify(simp)
    assert br.external.functional_derivative("src_Hdag", lab, "left").simplify(simp) == (-shd)
    assert br.external.functional_derivative("src_Hdag", lab, "right").simplify(simp) == shd
    pairs = {p.field: (p.sign, p.side) for p in br.pairs}
    assert pairs["Hdag"] == (1, "right") and pairs["Qbar"] == (-1, "left") and pairs["Q"] == (1, "left")
    assert pairs["H"] == (1, "left") and pairs["G"] == (1, "left") and pairs["gh_W"] == (1, "left")


def test_graded_leibniz_rule_of_s(sm, simp):
    """s(X Y) = (s X) Y + (-1)^{|X|} X (s Y) on X = Qbar (odd), Y = Q and on X = H (even), Y = Lbar."""
    s = BRST(sm).s
    for (nx, ny) in (("Qbar", "Q"), ("H", "Lbar"), ("gh_W", "W")):
        lx, ly = labels(sm, nx), labels(sm, ny)
        x, y = mono(sm, 1, Occurrence(nx, lx)), mono(sm, 1, Occurrence(ny, ly))
        sign = -1 if sm.fields[nx].odd else 1
        lhs = s(x * y)
        rhs = s(x) * y + (x * s(y)).scale(sign)
        assert lhs.simplify(simp) == rhs.simplify(simp), (nx, ny)


def test_vertex_function_leg_exchange(sm, simp):
    """Exchanging two odd legs flips V; moving an even leg does not."""
    S0 = BRST(sm).action
    legs = standard_legs(sm.fields, ("Qbar", "Q", "G"))
    v = S0.vertex_function(legs, simp)
    assert (S0.vertex_function([legs[1], legs[0], legs[2]], simp) + v).expand() == Expression.num(0)
    assert (S0.vertex_function([legs[0], legs[2], legs[1]], simp) - v).expand() == Expression.num(0)  # odd order kept
    assert (S0.vertex_function([legs[2], legs[1], legs[0]], simp) + v).expand() == Expression.num(0)
    legs = standard_legs(sm.fields, ("Hdag", "H", "B"))     # momenta: compare without eliminating one
    v = S0.vertex_function(legs, simp, momentum_conservation=False)
    w = S0.vertex_function([legs[2], legs[1], legs[0]], simp, momentum_conservation=False)
    assert (w - v).expand() == Expression.num(0) and v.expand() != Expression.num(0)


def test_conjugation(sm, simp):
    """(Qbar gamma5 Q)^dag = - Qbar gamma5 Q; (i Qbar gamma5 Q)^dag = + itself; (Qbar gammabar^mu Q B_mu)^dag = itself;
    (zeta_a f^{abc} omega^b omega^c)^dag = itself (zeta^dag = -zeta and the reversal of omega omega);
    (omegabar^a b^a)^dag = - itself."""
    qb, q = labels(sm, "Qbar"), labels(sm, "Q")
    d = tn.delta(qb[1], q[1]) * tn.delta(qb[2], q[2])
    g5 = mono(sm, tn.gamma5(qb[0], q[0]) * d, Occurrence("Qbar", qb), Occurrence("Q", q))
    assert g5.conjugate().simplify(simp) == (-g5).simplify(simp)
    ig5 = g5.scale(tn.I)
    assert ig5.conjugate().simplify(simp) == ig5.simplify(simp)
    mu = tn.fresh(tn.LORENTZ)
    cur = mono(sm, tn.gammabar(mu, qb[0], q[0]) * d, Occurrence("Qbar", qb), Occurrence("Q", q), Occurrence("B", (mu,)))
    assert cur.conjugate().simplify(simp) == cur.simplify(simp)
    a, b, c = labels(sm, "gh_G")[0], labels(sm, "gh_G")[0], labels(sm, "gh_G")[0]
    z = mono(sm, tn.f(a, b, c), Occurrence("src_gh_G", (a,)), Occurrence("gh_G", (b,)), Occurrence("gh_G", (c,)))
    assert z.conjugate() == z
    ab = labels(sm, "agh_G")[0]
    ob = mono(sm, 1, Occurrence("agh_G", (ab,)), Occurrence("b_G", (ab,)))
    assert ob.conjugate() == -ob


def test_index_spaces_are_never_promoted(sm, simp):
    """gbar . ghat = 0, gammabar . ghat = 0, gammahat . gbar = 0, ghat^mu_mu = d - 4, eps . ghat = 0;
    gamma5 anticommutes with gammabar and commutes with gammahat; the full metric splits into bar + hat."""
    mu, nu, rho, sg, tau = (tn.fresh(tn.LORENTZ) for _ in range(5))
    s, t, u = (tn.fresh(tn.SPINOR) for _ in range(3))
    z = Expression.num(0)
    assert simp(tn.gbar(mu, nu) * tn.ghat(nu, rho)) == z
    assert simp(tn.gammabar(mu, s, t) * tn.ghat(mu, nu)) == z
    assert simp(tn.gammahat(mu, s, t) * tn.gbar(mu, nu)) == z
    assert simp(tn.ghat(mu, mu)) == simp(tn.d - 4)
    assert simp(tn.eps(mu, nu, rho, sg) * tn.ghat(sg, tau)) == z
    assert simp(tn.gamma5(s, t) * tn.gammabar(mu, t, u) + tn.gammabar(mu, s, t) * tn.gamma5(t, u)) == z
    assert simp(tn.gamma5(s, t) * tn.gammahat(mu, t, u) - tn.gammahat(mu, s, t) * tn.gamma5(t, u)) == z
    assert simp(tn.g(mu, nu)) == simp(tn.gbar(mu, nu) + tn.ghat(mu, nu))
    assert simp(tn.p(1, mu) * tn.gbar(mu, nu) * tn.p(1, nu)) == tn.spbar(1, 1)
    assert simp(tn.p(2, mu) * tn.ghat(mu, nu) * tn.p(2, nu)) == tn.sphat(2, 2)
    assert simp(tn.p(1, mu) * tn.p(1, mu)) == simp(tn.spbar(1, 1) + tn.sphat(1, 1))
    assert tn.is_evanescent(simp(tn.gammahat(mu, s, t) * tn.p(1, mu)))
    assert not tn.is_evanescent(simp(tn.gammabar(mu, s, t) * tn.p(1, mu)))

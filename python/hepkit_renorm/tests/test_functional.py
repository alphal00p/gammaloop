"""Grassmann bookkeeping, vertex functions and the operations of :mod:`hepkit_renorm.functional`.

The reference signs are the report's (03_setup eq:V-src-gauge, eq:V-src-matter, 02_method sec:m-objects,
03b sec:st-signs) and those of the scout prototype ``e06_vertex_from_monomial.py``."""

import pytest
from symbolica import Expression

from hepkit_renorm import Density, Derivation, FieldImage, LocalFunctional, Occurrence, standard_legs
from hepkit_renorm import tensors as tn

I = tn.I


@pytest.fixture(scope="module")
def ctx(sm):
    return {"m": sm, "F": sm.fields, "simp": sm.simplifier(), "G": tn.factor_symbol("G"),
            "gs": sm.coupling("G"), "g": sm.coupling("W")}


def labels(ctx, name):
    return ctx["m"].slot_labels(name)


def rule(ctx, Y: LocalFunctional, names):
    """Feynman rule i V on the standard legs, simplified."""
    legs = standard_legs(ctx["F"], names)
    return legs, ctx["simp"](I * Y.vertex_function(legs, ctx["simp"]))


def same(ctx, a, b) -> bool:
    return ctx["simp"](a - b).expand() == Expression.num(0)


def test_canonical_order_carries_the_grassmann_sign(ctx):
    F = ctx["F"]
    q, qb, h = labels(ctx, "Q"), labels(ctx, "Qbar"), labels(ctx, "H")
    a = Density.monomial(F, 1, Occurrence("Qbar", qb), Occurrence("Q", q))
    b = Density.monomial(F, 1, Occurrence("Q", q), Occurrence("Qbar", qb))
    assert a == -b and not a.is_zero()
    c = Density.monomial(F, 1, Occurrence("H", h), Occurrence("Q", q))
    assert c == Density.monomial(F, 1, Occurrence("Q", q), Occurrence("H", h))


def test_grassmann_antisymmetry_kills_symmetric_ghost_pairs(ctx):
    F, G = ctx["F"], ctx["G"]
    a, b, c = tn.adj(G, 1), tn.adj(G, 2), tn.adj(G, 3)
    sym = Density.monomial(F, tn.delta(a, b), Occurrence("gh_G", (a,)), Occurrence("gh_G", (b,)))
    assert sym.is_zero()
    anti = Density.monomial(F, tn.f(c, a, b), Occurrence("gh_G", (a,)), Occurrence("gh_G", (b,)))
    assert len(anti) == 1 and anti.free_labels() == [c]


def test_equal_monomials_with_different_dummies_merge(ctx):
    F = ctx["F"]
    x, y = labels(ctx, "H"), labels(ctx, "H")
    two = Density.monomial(F, tn.delta(x[0], y[0]), Occurrence("Hdag", x), Occurrence("H", y))
    three = Density.monomial(F, tn.delta(x[0], y[0]), Occurrence("H", y), Occurrence("Hdag", x))
    assert len(two + three) == 1
    assert two.simplify(ctx["simp"]) == Density.monomial(F, 1, Occurrence("Hdag", x), Occurrence("H", x))


def test_source_vertex_signs(ctx):
    """Report eq:V-src-gauge and eq:V-src-matter (Feynman rules of the terms of eq:Lext)."""
    F, m, gs, g = ctx["F"], ctx["m"], ctx["gs"], ctx["g"]
    a, b, c = (tn.fresh(m.factors["G"].adjoint_space()) for _ in range(3))
    # -1/2 zeta_a fg^{abc} omega^b omega^c  ->  -i fg^{abc}
    Y = LocalFunctional.integrate(Density.monomial(
        F, -gs * tn.f(a, b, c) / 2, Occurrence("src_gh_G", (a,)), Occurrence("gh_G", (b,)), Occurrence("gh_G", (c,))))
    legs, r = rule(ctx, Y, ["src_gh_G", "gh_G", "gh_G"])
    (a1,), (a2,), (a3,) = (l.labels for l in legs)
    assert same(ctx, r, -I * gs * tn.f(a1, a2, a3))
    # rho^mu_a d_mu omega^a -> delta^{ab} p_mu with p the incoming ghost momentum, i.e. Gamma_{rho omega}(-q, q)
    # = -i q_mu = -Gamma_{omega rho}: on legs (rho(p1), omega(p2 = -p1)) the rule is -delta p1; on legs
    # (omega(p1), rho) the odd legs are exchanged and the rule is -delta p1 again
    mu = tn.fresh(tn.LORENTZ)
    Y = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("src_G", (mu, a)), Occurrence("gh_G", (a,), (mu,))))
    legs, r = rule(ctx, Y, ["src_G", "gh_G"])
    (m1, a1), (a2,) = (l.labels for l in legs)
    assert same(ctx, r, -tn.delta(a1, a2) * tn.p(1, m1))
    legs, r2 = rule(ctx, Y, ["gh_G", "src_G"])
    (a2,), (m1, a1) = (l.labels for l in legs)
    assert same(ctx, r2, -tn.delta(a1, a2) * tn.p(1, m1))
    # rho^mu_a fg^{abc} A^b_mu omega^c -> i fg^{abc} g_mu^nu
    nu = tn.fresh(tn.LORENTZ)
    Y = LocalFunctional.integrate(Density.monomial(
        F, gs * tn.f(a, b, c) * tn.g(mu, nu), Occurrence("src_G", (mu, a)), Occurrence("G", (nu, b)),
        Occurrence("gh_G", (c,))))
    legs, r = rule(ctx, Y, ["src_G", "G", "gh_G"])
    (m1, a1), (m2, a2), (a3,) = (l.labels for l in legs)
    assert same(ctx, r, I * gs * tn.f(a1, a2, a3) * tn.g(m1, m2))
    # i Rbar_f omega^a Tg^a P_c Psi_f -> -Tg^a P_c ; i omega^a Psibar_f P_cbar Tg^a R_f -> -P_cbar Tg^a
    s, t = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
    i, j, w = tn.fresh(m.fields["Q"].slots[1].space), tn.fresh(m.fields["Q"].slots[1].space), \
        tn.fresh(m.fields["Q"].slots[2].space)
    Y = LocalFunctional.integrate(Density.monomial(
        F, I * gs * tn.T(a, i, j) * tn.PL(s, t), Occurrence("src_Q", (s, i, w)), Occurrence("gh_G", (a,)),
        Occurrence("Q", (t, j, w))))
    legs, r = rule(ctx, Y, ["src_Q", "gh_G", "Q"])
    (s1, i1, w1), (a2,), (s3, i3, w3) = (l.labels for l in legs)
    ref = -gs * tn.T(a2, i1, i3) * tn.PL(s1, s3) * tn.delta(w1, w3)
    assert same(ctx, r, ref)
    # src_Q is even: moving it past the ghost changes nothing; exchanging the two odd legs flips the sign
    assert same(ctx, rule(ctx, Y, ["gh_G", "src_Q", "Q"])[1],
                _relabel_legs(ctx, ref, ["src_Q", "gh_G", "Q"], ["gh_G", "src_Q", "Q"]))
    flipped = rule(ctx, Y, ["Q", "gh_G", "src_Q"])[1]
    assert same(ctx, flipped, -_relabel_legs(ctx, ref, ["src_Q", "gh_G", "Q"], ["Q", "gh_G", "src_Q"]))
    Y = LocalFunctional.integrate(Density.monomial(
        F, I * gs * tn.PR(s, t) * tn.T(a, i, j), Occurrence("gh_G", (a,)), Occurrence("Qbar", (s, i, w)),
        Occurrence("src_Qbar", (t, j, w))))
    legs, r = rule(ctx, Y, ["gh_G", "Qbar", "src_Qbar"])
    (a1,), (s2, i2, w2), (s3, i3, w3) = (l.labels for l in legs)
    assert same(ctx, r, -tn.PR(s2, s3) * gs * tn.T(a1, i2, i3) * tn.delta(w2, w3))
    # i sigma^dag omega^a Tg^a_H H -> -Tg^a_H ;  -i omega^a H^dag Tg^a_H sigma -> +Tg^a_H   (sigma^dag = src_H)
    aw = tn.fresh(m.factors["W"].adjoint_space())
    x, y = tn.fresh(m.fields["H"].slots[0].space), tn.fresh(m.fields["H"].slots[0].space)
    Y = LocalFunctional.integrate(Density.monomial(
        F, I * g * tn.T(aw, x, y), Occurrence("src_H", (x,)), Occurrence("gh_W", (aw,)), Occurrence("H", (y,))))
    legs, r = rule(ctx, Y, ["src_H", "gh_W", "H"])
    (x1,), (a2,), (y3,) = (l.labels for l in legs)
    assert same(ctx, r, -g * tn.T(a2, x1, y3))
    Y = LocalFunctional.integrate(Density.monomial(
        F, -I * g * tn.T(aw, x, y), Occurrence("gh_W", (aw,)), Occurrence("Hdag", (x,)), Occurrence("src_Hdag", (y,))))
    legs, r = rule(ctx, Y, ["gh_W", "Hdag", "src_Hdag"])
    (a1,), (x2,), (y3,) = (l.labels for l in legs)
    assert same(ctx, r, g * tn.T(a1, x2, y3))


def _relabel_legs(ctx, expr, old_order, new_order):
    """Express a momentum-free leg-labelled expression in the labels of a permuted leg order (distinct fields)."""
    F = ctx["F"]
    old, new = standard_legs(F, old_order), standard_legs(F, new_order)
    by_field = {l.field: l for l in new}
    mapping = []
    for l in old:
        mapping += list(zip(l.labels, by_field[l.field].labels))
    return tn.rename(expr, mapping)


def test_functional_derivatives(ctx):
    F = ctx["F"]
    q, h = labels(ctx, "Q"), labels(ctx, "H")
    Fq = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("Qbar", q), Occurrence("Q", q)))
    L = labels(ctx, "Q")
    assert Fq.functional_derivative("Qbar", L) == Density.monomial(F, 1, Occurrence("Q", L))
    assert Fq.functional_derivative("Q", L) == Density.monomial(F, -1, Occurrence("Qbar", L))
    assert Fq.functional_derivative("Q", L, side="right") == Density.monomial(F, 1, Occurrence("Qbar", L))
    assert Fq.functional_derivative("Qbar", L, side="right") == Density.monomial(F, -1, Occurrence("Q", L))
    # Euler-Lagrange: d/dH^dag of int d H^dag . d H = - box H
    mu = tn.fresh(tn.LORENTZ)
    Fh = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("Hdag", h, (mu,)), Occurrence("H", h, (mu,))))
    Lh = labels(ctx, "H")
    nu = tn.fresh(tn.LORENTZ)
    assert Fh.functional_derivative("Hdag", Lh) == Density.monomial(F, -1, Occurrence("H", Lh, (nu, nu)))
    # report 03b sec:st-signs: for an even F and an odd field, dR F = - dL F
    m, aw = ctx["m"], tn.fresh(ctx["m"].factors["W"].adjoint_space())
    x, y = labels(ctx, "H")[0], labels(ctx, "H")[0]
    Fs = LocalFunctional.integrate(Density.monomial(
        F, -I * ctx["g"] * tn.T(aw, x, y), Occurrence("gh_W", (aw,)), Occurrence("Hdag", (x,)),
        Occurrence("src_Hdag", (y,))))
    Ls = labels(ctx, "src_Hdag")
    assert Fs.functional_derivative("src_Hdag", Ls, "right") == -Fs.functional_derivative("src_Hdag", Ls, "left")
    assert m.fields["src_Hdag"].is_source and m.fields["src_Hdag"].odd


def test_graded_derivation_is_a_left_derivation(ctx):
    F, m, gs = ctx["F"], ctx["m"], ctx["gs"]
    qlab = m.slot_labels("Q")
    a, j = tn.fresh(m.factors["G"].adjoint_space()), tn.fresh(F["Q"].slots[1].space)
    img_q = Density.monomial(F, I * gs * tn.T(a, qlab[1], j), Occurrence("gh_G", (a,)),
                             Occurrence("Q", (qlab[0], j, qlab[2])))
    qblab = m.slot_labels("Qbar")
    b, k = tn.fresh(m.factors["G"].adjoint_space()), tn.fresh(F["Q"].slots[1].space)
    img_qb = Density.monomial(F, I * gs * tn.T(b, k, qblab[1]), Occurrence("Qbar", (qblab[0], k, qblab[2])),
                              Occurrence("gh_G", (b,)))
    s = Derivation(1, {"Q": FieldImage(qlab, img_q), "Qbar": FieldImage(qblab, img_qb)})
    L = m.slot_labels("Q")
    A = Density.monomial(F, 1, Occurrence("Qbar", L))
    B = Density.monomial(F, 1, Occurrence("Q", L))
    assert s(A * B) == s(A) * B - A * s(B)  # (-1)^{|A|} with |A| = 1
    mu = tn.fresh(tn.LORENTZ)
    dB = Density.monomial(F, 1, Occurrence("Q", L, (mu,)))
    assert s(dB) == s(B).total_derivative(mu)  # s commutes with partial derivatives
    with pytest.raises(ValueError):
        Derivation(0, {"Q": FieldImage(qlab, img_q)})(B)  # image of the wrong parity


def test_integration_by_parts_and_total_derivatives(ctx):
    F, simp = ctx["F"], ctx["simp"]
    m = ctx["m"]
    mu = tn.fresh(tn.LORENTZ)
    a = tn.fresh(m.factors["G"].adjoint_space())
    rho_omega = Density.monomial(F, 1, Occurrence("src_G", (mu, a)), Occurrence("gh_G", (a,)))
    assert LocalFunctional.integrate(rho_omega.total_derivative(mu)).vanishes(simp)
    nonzero = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("src_G", (mu, a)),
                                                         Occurrence("gh_G", (a,), (mu,))))
    assert not nonzero.vanishes(simp)
    h = labels(ctx, "H")
    nu = tn.fresh(tn.LORENTZ)
    kin = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("Hdag", h, (nu,)), Occurrence("H", h, (nu,))))
    moved = kin.integrate_by_parts("Hdag")
    assert all(not o.derivatives for mono in moved for o in mono.occurrences if o.field == "Hdag")
    assert LocalFunctional.integrate(kin - moved).vanishes(simp)


def test_conjugation_is_an_involution_and_detects_antihermitian_terms(ctx):
    F, simp = ctx["F"], ctx["simp"]
    q = labels(ctx, "Q")
    mu = tn.fresh(tn.LORENTZ)
    s, t = q[0], tn.fresh(tn.SPINOR)
    kin = LocalFunctional.integrate(Density.monomial(
        F, I * tn.gamma(mu, s, t), Occurrence("Qbar", q), Occurrence("Q", (t,) + q[1:], (mu,))))
    assert kin.conjugate().conjugate() == kin
    assert LocalFunctional.integrate(kin - kin.conjugate()).vanishes(simp)  # hermitian up to IBP
    anti = LocalFunctional.integrate(kin.scale(I))
    assert not LocalFunctional.integrate(anti - anti.conjugate()).vanishes(simp)


def test_bmhv_projections(ctx):
    F, simp = ctx["F"], ctx["simp"]
    h = labels(ctx, "H")
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    trace = LocalFunctional.integrate(Density.monomial(F, tn.ghat(mu, nu) * tn.g(mu, nu),
                                                       Occurrence("Hdag", h), Occurrence("H", h)))
    # the closed hatted trace is d - 4: physical before d -> 4, zero in the four-dimensional limit
    assert trace.physical(simp) == trace.simplify(simp) and trace.evanescent(simp).is_zero()
    assert trace.four_dimensional_limit(simp).is_zero()
    kin = LocalFunctional.integrate(Density.monomial(F, 1, Occurrence("Hdag", h, (mu,)), Occurrence("H", h, (mu,))))
    a, b = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    e_h = LocalFunctional.integrate(Density.monomial(F, tn.ghat(a, b), Occurrence("Hdag", h, (a,)),
                                                     Occurrence("H", h, (b,))))  # report eq:evanescent E_H
    assert kin.evanescent(simp) == e_h.split_lorentz(simp)
    assert kin.physical(simp) + kin.evanescent(simp) == kin.split_lorentz(simp)
    assert kin.four_dimensional_limit(simp) == kin.physical(simp)

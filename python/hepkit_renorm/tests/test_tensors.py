from symbolica import Expression, S

from hepkit_renorm import tensors as tn

simp = tn.Simplifier({"G": (3, 8)}, n_flavours=1)
keep = tn.Simplifier({"G": (3, 8)}, expand_projectors=False)
mu, nu, rho = tn.lor(1), tn.lor(2), tn.lor(3)
s1, s2, s3, s4 = tn.spin(1), tn.spin(2), tn.spin(3), tn.spin(4)
X = S("hepkit_renorm_test::X", is_real=True)  # an external vector carrying a Lorentz label


def bar_gamma(m, a, b):
    n = tn.fresh(tn.LORENTZ)
    return tn.gamma(n, a, b) * tn.gbar(n, m)  # projector form, simplifies to gammabar(m, a, b)


def hat_gamma(m, a, b):
    return tn.gammahat(m, a, b)


def test_metric_algebra():
    assert simp(tn.g(mu, mu)) == tn.d
    assert simp(tn.gbar(mu, nu) * tn.gbar(nu, mu)) == Expression.num(4)
    assert simp(tn.ghat(mu, nu) * tn.ghat(nu, mu)) == tn.d - 4
    assert simp(tn.gbar(mu, nu) * tn.ghat(nu, rho)) == Expression.num(0)
    assert simp(tn.g(mu, nu) * tn.gbar(nu, rho) * X(mu)) == simp(tn.gbar(mu, rho) * X(mu))
    assert simp(tn.eps(mu, nu, rho, tn.lor(4)) * tn.ghat(tn.lor(4), tn.lor(5))) == Expression.num(0)
    assert simp(tn.p(1, mu) * tn.p(2, mu)) == tn.spbar(1, 2) + tn.sphat(1, 2)  # canonical mode splits p.p
    assert keep(tn.p(1, mu) * tn.p(2, mu)) == tn.sp(1, 2)
    assert simp(tn.p(1, mu) * tn.gbar(mu, nu) * tn.p(1, nu)) == tn.spbar(1, 1)
    assert simp(tn.gamma(nu, s1, s2) * tn.gbar(nu, mu)) == tn.gammabar(mu, s1, s2)
    assert simp(tn.gamma(mu, s1, s2)) == tn.gammabar(mu, s1, s2) + tn.gammahat(mu, s1, s2)
    assert simp(tn.gammahat(mu, s1, s2) * tn.gbar(mu, nu)) == 0
    assert simp(tn.delta(tn.spin(1), tn.spin(1))) == Expression.num(4)


def test_gamma5_with_barred_and_hatted_gammas():
    # {gbar^mu, gamma5} = 0, [ghat^mu, gamma5] = 0, gamma5^2 = 1 (report eq:g5comm)
    assert simp(tn.gamma5(s1, s2) * bar_gamma(mu, s2, s3) + bar_gamma(mu, s1, s2) * tn.gamma5(s2, s3)) == 0
    assert simp(tn.gamma5(s1, s2) * hat_gamma(mu, s2, s3) - hat_gamma(mu, s1, s2) * tn.gamma5(s2, s3)) == 0
    assert simp(tn.gamma5(s1, s2) * tn.gamma5(s2, s3)) == tn.delta(s1, s3)
    # {gamma^mu, gamma5} = 2 ghat^mu_nu gamma^nu gamma5
    lhs = tn.gamma5(s1, s2) * tn.gamma(mu, s2, s3) + tn.gamma(mu, s1, s2) * tn.gamma5(s2, s3)
    assert simp(lhs - 2 * hat_gamma(mu, s1, s2) * tn.gamma5(s2, s3)) == 0


def test_sandwich_rules():
    # PR gamma^mu PL = gbar^mu PL and PL gamma^mu PL = ghat^mu PL (report eq:sandwich), in both modes
    for sim in (simp, keep):
        assert sim(tn.PR(s1, s2) * tn.gamma(mu, s2, s3) * tn.PL(s3, s4) - bar_gamma(mu, s1, s2) * tn.PL(s2, s4)) == 0
        assert sim(tn.PL(s1, s2) * tn.gamma(mu, s2, s3) * tn.PL(s3, s4) - hat_gamma(mu, s1, s2) * tn.PL(s2, s4)) == 0
        assert sim(tn.PL(s1, s2) * tn.PR(s2, s3)) == 0


def test_conjugation():
    g_real = S("hepkit_renorm_test::gr", is_real=True)
    y = S("hepkit_renorm_test::y")
    c = (Expression.num(1) + tn.I) * g_real * y * tn.gamma(mu, s1, s2) * tn.gamma5(s2, s3) * tn.PL(s3, s4)
    conj = tn.conjugate(c)
    expected = ((Expression.num(1) - tn.I) * g_real * y.conj() * tn.gamma(mu, s2, s1) * (-tn.gamma5(s3, s2))
                * tn.PR(s4, s3))
    assert (conj - expected).expand() == 0
    assert (tn.conjugate(conj) - c).expand() == 0
    a, i, j = tn.adj(tn.factor_symbol("G"), 1), tn.fund(tn.factor_symbol("G"), 1), tn.fund(tn.factor_symbol("G"), 2)
    assert tn.conjugate(tn.T(a, i, j)) == tn.T(a, j, i)


def test_canonical_relabelling_and_labels():
    a = tn.g(tn.lor(5), tn.lor(9)) * X(tn.lor(5)) * X(tn.lor(9)) * tn.gamma(tn.lor(7), s1, s2) * X(tn.lor(7))
    b = tn.g(tn.lor(11), tn.lor(3)) * X(tn.lor(3)) * X(tn.lor(11)) * tn.gamma(tn.lor(8), s1, s2) * X(tn.lor(8))
    assert tn.canonical_relabel(a) == tn.canonical_relabel(b)
    assert {x.to_canonical_string() for x in tn.free_labels(a)} == {s1.to_canonical_string(), s2.to_canonical_string()}
    assert tn.label_space(tn.adj(tn.factor_symbol("W"), 3)) == "adjoint:W"
    assert tn.contains_head(a, tn.gamma) and not tn.contains_head(a, tn.ghat)
    f1 = tn.fresh(tn.LORENTZ)
    assert f1 != tn.fresh(tn.LORENTZ)

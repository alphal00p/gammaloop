"""BRST transformations, gauge fermion, sources, the Slavnov-Taylor operator, b_d / b and Delta_1 of the one-family
Standard Model (gate V2 of the validation ladder).

References: report 03_setup sec:brst (eq:brst, eq:brst-matter, eq:Lgf, eq:Lext, eq:gf-exact, eq:ST, eq:Phistar,
eq:b-sources, eq:linear-ids, eq:rhotilde), sec:delta (eq:Delta1, eq:V-delta), sec:feynman-rules (eq:V-src-gauge,
eq:V-src-matter), 03b sec:st-signs, sec:dictionary (examples 1-3), sec:btable (eq:btable-gauge, eq:btable-fermion),
02_method eq:m-S-exact, and the plane-wave engine tests ``SCRATCHPAD/st/test_signs.py``, ``test_b2.py``,
``test_b2_control.py`` and ``test_blin.py``, which are reproduced here on the generated objects.

Zero tests are exact and structural: vertex functions after the BMHV simplifier, with the group tensors evaluated
over a prime field where Lie-algebra identities are needed (:func:`hepkit_renorm.brst.vanishes`). A density with
free labels is tested through its pairing with a field (:func:`hepkit_renorm.brst.paired`).
"""

from fractions import Fraction

import pytest
from symbolica import Expression, S

from hepkit_renorm import Density, Derivation, FieldImage, LocalFunctional, Occurrence, standard_legs
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST, SourcePair, paired, vanishes
from hepkit_renorm.groups import generator_expr

I = tn.I
_w = S("w_")
FERMIONS = ("Q", "U", "D", "L", "E")
PAIRS_OF_FACTORS = [(f, r) for f, rs in {"Q": "GWB", "U": "GB", "D": "GB", "L": "WB", "E": "B"}.items() for r in rs]


@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def simp(sm):
    return sm.simplifier()


def num(q) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


def zero(sm, functional: LocalFunctional) -> bool:
    return vanishes(sm, LocalFunctional.integrate(functional))


def equal(sm, a: Density, b: Density) -> bool:
    return zero(sm, a - b)


def zero4(sm, functional: Density) -> bool:
    """Zero in four dimensions: the four-dimensional limit (P, d -> 4) of the functional vanishes."""
    return zero(sm, functional.four_dimensional_limit(sm.simplifier()))


_FINGERPRINTERS: dict = {}


def nonzero(sm, functional: Density, dimension: int = 6) -> bool:
    """Exact certificate that a functional is NOT zero: a non-vanishing plane-wave fingerprint (``zero`` returning
    False is not one, since the simplifier does not know every identity). Used by the negative controls."""
    from hepkit_renorm.fingerprint import Fingerprinter

    key = (id(sm), dimension)
    if key not in _FINGERPRINTERS:
        _FINGERPRINTERS[key] = Fingerprinter(sm, dimension=dimension)
    return not _FINGERPRINTERS[key].vanishes(LocalFunctional.integrate(functional))


def nonzero4(sm, functional: Density) -> bool:
    """Non-zero in four dimensions (exact certificate, see :func:`nonzero`)."""
    return nonzero(sm, functional.four_dimensional_limit(sm.simplifier()), dimension=4)


def same_expr(sm, value, reference, last_momentum=None) -> bool:
    simp = sm.simplifier()
    if last_momentum is not None:
        others = sum((tn.p(k, _w) for k in range(1, last_momentum)), Expression.num(0))
        reference = reference.replace(tn.p(last_momentum, _w), -others)
    return simp(value - reference).expand() == Expression.num(0)


def quantum_numbers(sm, density: Density) -> set[tuple[int, Fraction]]:
    """{(ghost number, mass dimension)} of the monomials (derivatives count one)."""
    out = set()
    for mono in density.monomials():
        gh = sum(sm.fields[o.field].ghost_number for o in mono.occurrences)
        dim = sum(sm.fields[o.field].dimension + len(o.derivatives) for o in mono.occurrences)
        out.add((gh, dim))
    return out


def generator(sm, multiplet: str, factor: str, a, i_labels, j_labels):
    """Coupling-weighted generator cal(T)^a_f between the group labels of two legs (deltas on the other slots)."""
    fac = sm.factors[factor]
    out = sm.coupling(factor)
    slots = sm.fields[multiplet].slots
    group = [k for k, sl in enumerate(slots) if sl.space.kind in ("adjoint", "fund")]
    for pos, k in enumerate(group):
        i, j = i_labels[pos], j_labels[pos]
        out = out * (generator_expr(fac, slots[k].rep, a, i, j) if slots[k].factor == factor else tn.delta(i, j))
    if fac.abelian:
        out = out * num(sm.fields[multiplet].charges.get(factor, 0))
    return out


def projector(sm, f, physical=True):
    left = sm.fermions[f].chirality == "L"
    return tn.PL if left == physical else tn.PR


# ---------------------------------------------------------------------------------------------- s on fields
def test_brst_transformations_of_the_fields(sm, brst, simp):
    """report eq:brst and eq:brst-matter: s A = D omega, s omega = -1/2 fg omega omega, s omegabar = b, s b = 0,
    s Psi = i omega cal(T) P_c Psi, s Psibar = i Psibar P_cbar cal(T) omega, s H = i omega cal(T) H,
    s H^dag = -i H^dag cal(T) omega; sources are invariant."""
    F, gs, g, gp = sm.fields, sm.coupling("G"), sm.coupling("W"), sm.coupling("B")

    def mono(c, *occ):
        return Density.monomial(F, c, *occ)

    # gauge field and ghost of SU(2)
    mu, a = tn.fresh(tn.LORENTZ), tn.fresh(sm.factors["W"].adjoint_space())
    b, c = tn.fresh(sm.factors["W"].adjoint_space()), tn.fresh(sm.factors["W"].adjoint_space())
    ref = mono(1, Occurrence("gh_W", (a,), (mu,))) + mono(g * tn.f(a, b, c), Occurrence("W", (mu, b)),
                                                          Occurrence("gh_W", (c,)))
    assert brst.variation("W", (mu, a)).simplify(simp) == ref.simplify(simp)
    ref = mono(-g * tn.f(a, b, c) / 2, Occurrence("gh_W", (b,)), Occurrence("gh_W", (c,)))
    assert brst.variation("gh_W", (a,)) == ref
    assert brst.variation("B", (mu,)) == mono(1, Occurrence("gh_B", (), (mu,)))
    assert brst.variation("gh_B", ()).is_zero()
    assert brst.variation("agh_W", (a,)) == mono(1, Occurrence("b_W", (a,)))
    assert brst.variation("b_W", (a,)).is_zero()
    assert all(brst.variation(n, sm.slot_labels(n)).is_zero() for n, fi in F.items() if fi.is_source)
    # left-handed quark doublet: s Q = i omega^a cal(T)^a PL Q over G, W, B
    s1, c1, w1 = sm.slot_labels("Q")
    u, c2, w2 = tn.fresh(tn.SPINOR), tn.fresh(F["Q"].slots[1].space), tn.fresh(F["Q"].slots[2].space)
    aG, aW = tn.fresh(sm.factors["G"].adjoint_space()), tn.fresh(sm.factors["W"].adjoint_space())
    ref = (mono(I * gs * tn.T(aG, c1, c2) * tn.PL(s1, u), Occurrence("gh_G", (aG,)), Occurrence("Q", (u, c2, w1)))
           + mono(I * g * tn.T(aW, w1, w2) * tn.PL(s1, u), Occurrence("gh_W", (aW,)), Occurrence("Q", (u, c1, w2)))
           + mono(I * gp / 6 * tn.PL(s1, u), Occurrence("gh_B", ()), Occurrence("Q", (u, c1, w1))))
    assert brst.variation("Q", (s1, c1, w1)).simplify(simp) == ref.simplify(simp)
    # right-handed up quark, conjugate: s Ubar = i Ubar PL cal(T) omega (ghost written to the right)
    s1, c1 = sm.slot_labels("Ubar")
    ref = (mono(I * gs * tn.PL(u, s1) * tn.T(aG, c2, c1), Occurrence("Ubar", (u, c2)), Occurrence("gh_G", (aG,)))
           + mono(I * gp * 2 / 3 * tn.PL(u, s1), Occurrence("Ubar", (u, c1)), Occurrence("gh_B", ())))
    assert brst.variation("Ubar", (s1, c1)).simplify(simp) == ref.simplify(simp)
    # Higgs doublet and its conjugate
    (i1,), j = sm.slot_labels("H"), tn.fresh(F["H"].slots[0].space)
    ref = (mono(I * g * tn.T(aW, i1, j), Occurrence("gh_W", (aW,)), Occurrence("H", (j,)))
           + mono(I * gp / 2, Occurrence("gh_B", ()), Occurrence("H", (i1,))))
    assert brst.variation("H", (i1,)).simplify(simp) == ref.simplify(simp)
    ref = (mono(-I * g * tn.T(aW, j, i1), Occurrence("Hdag", (j,)), Occurrence("gh_W", (aW,)))
           + mono(-I * gp / 2, Occurrence("Hdag", (i1,)), Occurrence("gh_B", ())))
    assert brst.variation("Hdag", (i1,)).simplify(simp) == ref.simplify(simp)


def test_only_the_physical_chirality_transforms(sm, brst, simp):
    """s chi_f = 0: the image of P_cbar Psi_f (and of Psibar_f P_c) vanishes (report sec:brst)."""
    for f in FERMIONS:
        for name, physical in ((f, False), (f + "bar", True)):
            labels = sm.slot_labels(name)
            u = tn.fresh(tn.SPINOR)
            image = brst.variation(name, (u,) + labels[1:])
            head = projector(sm, f, physical)
            proj = head(labels[0], u) if name == f else head(u, labels[0])
            assert image.scale(proj).simplify(simp).is_zero(), name
            assert not image.scale(projector(sm, f, not physical)(labels[0], u) if name == f else
                                   projector(sm, f, not physical)(u, labels[0])).simplify(simp).is_zero()


def test_s_is_nilpotent_on_every_field(sm, brst):
    """s^2 = 0 on every field (report sec:brst), up to the Lie algebra, tested through int X^* s^2 X."""
    for name, info in sm.fields.items():
        if info.is_source:
            continue
        labels = sm.slot_labels(name)
        s2 = brst.s(brst.variation(name, labels))
        if info.source is None:
            assert s2.is_zero(), name
        else:
            assert zero(sm, paired(sm, info.source, labels, s2)), name
    # the control: s^2 Q = 0 needs the Lie algebra; without the -1/2 fg omega omega of s omega it fails
    images = {k: v for k, v in brst.s.images.items() if k != "gh_G"}
    s_broken = Derivation(1, images)
    labels = sm.slot_labels("Q")
    assert nonzero(sm, paired(sm, "src_Q", labels, s_broken(s_broken(brst.variation("Q", labels)))))


# ---------------------------------------------------------------------------------------------- gauge fermion
def test_gauge_fixing_and_ghost_terms(sm, brst):
    """L_gf = s[omegabar (d.A + xi/2 b)] = b^a d.A^a + xi/2 b^a b^a + d omegabar^a (D omega)^a (report eq:Lgf)."""
    F, xi = sm.fields, sm.gauge_fixing.xi
    ref = Density(F)
    for r, fac in sm.factors.items():
        mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
        adj = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        ref = ref + Density.monomial(F, 1, Occurrence(f"b_{r}", adj), Occurrence(r, (mu,) + adj, (mu,)))
        ref = ref + Density.monomial(F, xi / 2, Occurrence(f"b_{r}", adj), Occurrence(f"b_{r}", adj))
        Dw = brst.variation(r, (nu,) + adj)  # (D_nu omega)^a
        ref = ref + Density.monomial(F, 1, Occurrence(f"agh_{r}", adj, (nu,))) * Dw
    assert equal(sm, brst.gauge_fixing, ref)


def test_gauge_fermion_reproduces_eq_gf_exact(sm, brst, simp):
    """The source part of the gauge fermion is - rho A + zeta omega + sum_f (Rbar Psi - Psibar R) - sigma^dag H
    + H^dag sigma (report eq:gf-exact); it is anti-hermitian, and L_ext, L_gf and S_0 are hermitian."""
    F = sm.fields
    ref = Density(F)
    for name, sign in [("G", -1), ("W", -1), ("B", -1), ("gh_G", 1), ("gh_W", 1), ("gh_B", 1), ("H", -1),
                       *[(f, 1) for f in FERMIONS]]:
        labels = sm.slot_labels(name)
        ref = ref + Density.monomial(F, sign, Occurrence(F[name].source, labels), Occurrence(name, labels))
    for name, sign in [("Hdag", 1), *[(f + "bar", -1) for f in FERMIONS]]:
        labels = sm.slot_labels(name)
        ref = ref + Density.monomial(F, sign, Occurrence(name, labels), Occurrence(F[name].source, labels))
    assert brst.gauge_fermion_sources.simplify(simp) == ref.simplify(simp)
    psi = brst.gauge_fermion
    assert zero(sm, psi + psi.conjugate())
    for part in (brst.external, brst.gauge_fixing, brst.action):
        assert zero(sm, part - part.conjugate())
    assert quantum_numbers(sm, psi) == {(-1, Fraction(3))}
    assert quantum_numbers(sm, brst.external) == {(0, Fraction(4))}
    assert {gh for gh, _ in quantum_numbers(sm, brst.action)} == {0}


def _src_rule(sm, brst, names):
    legs = standard_legs(sm.fields, names)
    return [leg.labels for leg in legs], sm.simplifier()(I * brst.external.vertex_function(legs, sm.simplifier()))


@pytest.mark.parametrize("r", ["G", "W", "B"])
def test_source_vertices_of_the_gauge_sector(sm, brst, r):
    """report eq:V-src-gauge from the derived L_ext: rho gh(p): delta p_mu; rho A gh: i fg g_mu^nu;
    zeta gh gh: -i fg (Feynman rules, momenta incoming)."""
    c, abelian = sm.coupling(r), sm.factors[r].abelian
    (m1, *a1), (*a2,) = _src_rule(sm, brst, [f"src_{r}", f"gh_{r}"])[0]
    delta = Expression.num(1) if abelian else tn.delta(a1[0], a2[0])
    assert same_expr(sm, _src_rule(sm, brst, [f"src_{r}", f"gh_{r}"])[1], delta * tn.p(2, m1), 2)
    if abelian:
        assert brst.external.vertex_function(standard_legs(sm.fields, ["src_gh_B", "gh_B", "gh_B"])) == 0
        return
    (m1, a1), (m2, a2), (a3,) = _src_rule(sm, brst, [f"src_{r}", r, f"gh_{r}"])[0]
    assert same_expr(sm, _src_rule(sm, brst, [f"src_{r}", r, f"gh_{r}"])[1], I * c * tn.f(a1, a2, a3) * tn.g(m1, m2))
    (a1,), (a2,), (a3,) = _src_rule(sm, brst, [f"src_gh_{r}", f"gh_{r}", f"gh_{r}"])[0]
    assert same_expr(sm, _src_rule(sm, brst, [f"src_gh_{r}", f"gh_{r}", f"gh_{r}"])[1], -I * c * tn.f(a1, a2, a3))


@pytest.mark.parametrize("f, r", PAIRS_OF_FACTORS)
def test_source_vertices_of_the_fermions(sm, brst, f, r):
    """report eq:V-src-matter: Rbar_f gh^a Psi_f: -cal(T)^a P_c ; gh^a Psibar_f R_f: -P_cbar cal(T)^a."""
    (s1, *g1), (*a2,), (s3, *g3) = _src_rule(sm, brst, [f"src_{f}", f"gh_{r}", f])[0]
    a = a2[0] if a2 else None
    ref = -generator(sm, f, r, a, g1, g3) * projector(sm, f)(s1, s3)
    assert same_expr(sm, _src_rule(sm, brst, [f"src_{f}", f"gh_{r}", f])[1], ref)
    (*a1,), (s2, *g2), (s3, *g3) = _src_rule(sm, brst, [f"gh_{r}", f"{f}bar", f"src_{f}bar"])[0]
    a = a1[0] if a1 else None
    ref = -projector(sm, f, physical=False)(s2, s3) * generator(sm, f, r, a, g2, g3)
    assert same_expr(sm, _src_rule(sm, brst, [f"gh_{r}", f"{f}bar", f"src_{f}bar"])[1], ref)


@pytest.mark.parametrize("r", ["W", "B"])
def test_source_vertices_of_the_higgs(sm, brst, r):
    """report eq:V-src-matter: sigma^dag gh^a H: -cal(T)^a_H ; gh^a H^dag sigma: +cal(T)^a_H."""
    (i1,), (*a2,), (j3,) = _src_rule(sm, brst, ["src_H", f"gh_{r}", "H"])[0]
    a = a2[0] if a2 else None
    assert same_expr(sm, _src_rule(sm, brst, ["src_H", f"gh_{r}", "H"])[1], -generator(sm, "H", r, a, [i1], [j3]))
    (*a1,), (i2,), (j3,) = _src_rule(sm, brst, [f"gh_{r}", "Hdag", "src_Hdag"])[0]
    a = a1[0] if a1 else None
    assert same_expr(sm, _src_rule(sm, brst, [f"gh_{r}", "Hdag", "src_Hdag"])[1], generator(sm, "H", r, a, [i2], [j3]))


# ---------------------------------------------------------------------------------------------- the ST operator
def test_source_pairs_reproduce_eq_ST(sm, brst):
    """The derived pairs are those of report eq:ST: left derivatives with a plus sign for A, omega, Psi, H, a minus
    sign for Psibar, the right derivative for H^dag; delta' S_0/delta Phi^* = s Phi (eq:ST-compact)."""
    expected = {}
    for r in sm.factors:
        expected[r] = expected[f"gh_{r}"] = (1, "left", 1)
    for f in FERMIONS:
        expected[f], expected[f + "bar"] = (1, "left", 1), (-1, "left", -1)
    expected["H"], expected["Hdag"] = (1, "left", 1), (1, "right", -1)
    assert {p.field: (p.sign, p.side, p.left_sign) for p in brst.pairs} == expected
    assert all(sm.fields[p.source].source_of == p.field for p in brst.pairs)
    assert brst.linear_fields == ("agh_G", "agh_W", "agh_B")
    simp = sm.simplifier()
    for p in brst.pairs:
        labels = sm.slot_labels(p.field)
        first = brst.action.functional_derivative(p.source, labels, p.side).scale(p.sign)
        assert first.simplify(simp) == brst.variation(p.field, labels).simplify(simp), p.field


def test_st_of_s0_is_s_s0_and_the_draft_operator_fails(sm, brst):
    """ST(S_0) = s S_0 (report eq:ST-S0; st/test_signs.py (a)). The draft forms of 03b sec:st-signs,
    (dR F/dPsi)(dF/dRbar) = -(dF/dRbar)(dL F/dPsi) and (dF/dH^dag)(dL F/dsigma) = (dL F/dsigma)(dF/dH^dag),
    give ST(S_0) != s S_0 on classes with fermion or H^dag legs."""
    st = brst.slavnov_taylor(brst.action)
    assert zero(sm, st - brst.breaking)
    assert zero4(sm, st) and not zero(sm, st)  # ST(S_0) = 0 in four dimensions only
    draft_fermion = [SourcePair(p.field, p.source, -p.sign, p.side, -p.left_sign) if p.field == "Q" else p
                     for p in brst.pairs]
    draft_higgs = [SourcePair(p.field, p.source, 1, "left", 1) if p.field == "Hdag" else p for p in brst.pairs]
    for pairs, legs in ((draft_fermion, ("Q", "Qbar")), (draft_higgs, ("Hdag",))):
        bad = LocalFunctional.integrate(brst.slavnov_taylor(brst.action, pairs) - brst.breaking)
        assert nonzero(sm, bad)
        assert all(any(x in c for x in legs) for c in bad.leg_classes() if not zero(sm, _class_part(bad, c)))
    with pytest.raises(ValueError):
        brst.slavnov_taylor(brst.gauge_fermion)


def _class_part(F: LocalFunctional, cls) -> LocalFunctional:
    return LocalFunctional.from_monomials(F.fields, [m for m in F.monomials() if m.fields() == cls])


def test_gauge_fermion_signs_matter(sm, brst):
    """report 03b sec:st-signs: with the opposite sign of Psibar_f R_f or of H^dag sigma in eq:gf-exact, the
    operator of eq:ST gives ST(S_0) != s S_0."""
    F = sm.fields
    for name in ("Qbar", "Hdag"):
        labels = sm.slot_labels(name)
        sign = -1 if name == "Qbar" else 1
        term = Density.monomial(F, sign, Occurrence(name, labels), Occurrence(F[name].source, labels))
        flipped = LocalFunctional.integrate(brst.gauge_fermion_sources - term.scale(2))
        action = LocalFunctional.integrate(brst.classical + brst.gauge_fixing + brst.s(flipped))
        assert nonzero(sm, brst.slavnov_taylor(action, brst.pairs) - brst.s(action)), name


def test_st_of_s0_plus_a_non_invariant_source_free_term(sm, brst):
    """ST(S_0 + T) - ST(S_0) = s T for a source-free, non-invariant T (st/test_signs.py (b))."""
    F = sm.fields
    s1, c1, w1 = sm.slot_labels("Q")
    t, mu, nu = tn.fresh(tn.SPINOR), tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    u = tn.fresh(tn.SPINOR)
    T = Density.monomial(F, num("3/10") * tn.gamma(mu, s1, u) * sm.coupling("B") / 6 * tn.PL(u, t),
                         Occurrence("Qbar", (s1, c1, w1)), Occurrence("Q", (t, c1, w1)), Occurrence("B", (mu,)))
    (i,) = sm.slot_labels("H")
    a = tn.fresh(sm.factors["W"].adjoint_space())
    for r, adj in (("W", (a,)), ("B", ())):
        T = T + Density.monomial(F, num("7/10") * tn.g(mu, nu), Occurrence("Hdag", (i,)), Occurrence("H", (i,)),
                                 Occurrence(r, (mu,) + adj), Occurrence(r, (nu,) + adj))
    T = LocalFunctional.integrate(T)
    assert sm.rigid_defects(T, _field()) == [] and not zero(sm, brst.s(T))  # rigid, but not BRST invariant
    st = brst.slavnov_taylor(LocalFunctional.integrate(brst.action + T))
    assert zero(sm, st - brst.breaking - brst.s(T))


def _field():
    from hepkit_renorm.exact import default_field

    return default_field()


def _source_test_functional(sm) -> LocalFunctional:
    """An even, non-invariant, non-hermitian G that contains every kind of source (st/test_blin.py)."""
    F = sm.fields
    aG = [tn.fresh(sm.factors["G"].adjoint_space()) for _ in range(3)]
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    G = Density.monomial(F, num("2/5") * tn.f(*aG), Occurrence("src_gh_G", (aG[0],)), Occurrence("gh_G", (aG[1],)),
                         Occurrence("gh_G", (aG[2],)))
    G = G + Density.monomial(F, num("9/10"), Occurrence("src_G", (mu, aG[0])), Occurrence("gh_G", (aG[0],), (mu,)))
    s, c, w = sm.slot_labels("Q")
    t, j = tn.fresh(tn.SPINOR), tn.fresh(F["Q"].slots[1].space)
    G = G + Density.monomial(F, num("3/5") * I * tn.T(aG[0], c, j) * tn.PL(s, t), Occurrence("src_Q", (s, c, w)),
                             Occurrence("gh_G", (aG[0],)), Occurrence("Q", (t, j, w)))
    s, c = sm.slot_labels("U")
    G = G + Density.monomial(F, I / 2 * tn.PL(s, t), Occurrence("gh_B", ()), Occurrence("Ubar", (s, c)),
                             Occurrence("src_Ubar", (t, c)))
    (i,), k = sm.slot_labels("H"), tn.fresh(F["H"].slots[0].space)
    aW = tn.fresh(sm.factors["W"].adjoint_space())
    G = G + Density.monomial(F, num("4/5") * I * tn.T(aW, i, k), Occurrence("src_H", (i,)),
                             Occurrence("gh_W", (aW,)), Occurrence("H", (k,)))
    G = G + Density.monomial(F, num("3/10") * I, Occurrence("gh_B", ()), Occurrence("Hdag", (i,)),
                             Occurrence("src_Hdag", (i,)))
    G = G + Density.monomial(F, num("7/10") * tn.g(mu, nu), Occurrence("Hdag", (i,)), Occurrence("H", (i,)),
                             Occurrence("B", (mu,)), Occurrence("B", (nu,)))
    return LocalFunctional.integrate(G)


def test_linearised_operator(sm, brst):
    """b_d G = ST(S_0 + G) - ST(S_0) - ST(G) for an even G containing every source (st/test_blin.py), and
    ST(S_0 + G) = Delta_1 + b_d G + S_q(G) (report eq:m-S-exact)."""
    G = _source_test_functional(sm)
    full = brst.slavnov_taylor(LocalFunctional.integrate(brst.action + G))
    bG = brst.b_d(G)
    assert nonzero(sm, bG)
    assert zero(sm, full - brst.slavnov_taylor(brst.action) - brst.slavnov_taylor(G) - bG)
    assert zero(sm, full - brst.breaking - bG - brst.quadratic_part(G))
    assert nonzero(sm, brst.quadratic_part(G))


def test_b_d_acts_as_s_on_fields_and_as_the_equations_of_motion_on_sources(sm, brst, simp):
    """report eq:b-sources: b_d Phi = s Phi (derived, not assumed) and b_d Phi^* = s_Phi dL S_0/dPhi."""
    for name, img in brst.s.images.items():
        bd = brst.b_d.images[name]
        assert bd.density.rename(list(zip(bd.labels, img.labels))).simplify(simp) == img.density.simplify(simp), name
    for p in brst.pairs:
        img = brst.b_d.images[p.source]
        eom = brst.action.functional_derivative(p.field, img.labels).scale(p.left_sign)
        assert img.density == eom, p.source
    assert {p.source: p.left_sign for p in brst.pairs if p.left_sign < 0} == {
        "src_Qbar": -1, "src_Ubar": -1, "src_Dbar": -1, "src_Lbar": -1, "src_Ebar": -1, "src_Hdag": -1}


def test_b_squared(sm, brst):
    """b^2 = 0 in four dimensions on every field and source (st/test_b2.py). In d dimensions b_d^2 vanishes on the
    fields and on rho, sigma, sigma^dag; on zeta, Rbar and R it is -s_Phi dL Delta_1/dPhi (Delta_1 depends on
    omega, Psi and Psibar only; report 03b sec:evanescent, 'closed directions in d dimensions')."""
    simp = sm.simplifier()
    for p in brst.pairs:
        for target, partner in ((p.field, p.source), (p.source, p.field)):
            img = brst.b.images.get(target)
            if img is None:
                assert brst.b_d.images.get(target) is None and target == "gh_B"
                continue
            probe = paired(sm, partner, img.labels, brst.b(img.density))
            assert zero(sm, probe.four_dimensional_limit(simp)), target
            bd = brst.b_d.images[target]
            d_dim = brst.b_d(bd.density)
            if target == p.source:
                d_dim = d_dim + brst.breaking.functional_derivative(p.field, bd.labels).scale(p.left_sign)
            assert zero(sm, paired(sm, partner, bd.labels, d_dim)), target
    labels = brst.b_d.images["src_Q"].labels
    assert nonzero(sm, paired(sm, "Q", labels, brst.b_d(brst.b_d.images["src_Q"].density)))


@pytest.mark.parametrize("source", ["src_Hdag", "src_Qbar", "src_L", "src_Q", "src_H", "src_W", "src_gh_G"])
def test_b_squared_detects_a_sign_flip(sm, brst, source):
    """Flipping the sign of b on one source gives b^2 != 0 on that source (st/test_b2_control.py)."""
    simp = sm.simplifier()
    images = dict(brst.b.images)
    img = images[source]
    images[source] = FieldImage(img.labels, -img.density)
    b = Derivation(1, images)
    probe = paired(sm, sm.fields[source].source_of, img.labels, b(images[source].density))
    assert nonzero(sm, probe.four_dimensional_limit(simp))


# ---------------------------------------------------------------------------------------------- Delta_1
def test_breaking_is_the_evanescent_fermion_kinetic_variation(sm, brst):
    """Delta_1 = sum_f int omega^a [Psibar cal(T)^a ghat^mu P_cbar d_mu Psi + (d_mu Psibar) ghat^mu P_c cal(T)^a Psi]
    (report eq:Delta1, with its sign); P Delta_1 = 0; ghost number 1, dimension 5; rigidly invariant."""
    F = sm.fields
    ref = Density(F)
    for f in FERMIONS:
        info = F[f]
        for r in sm.factors:
            if r not in sm.fermions[f].reps:
                continue
            fac = sm.factors[r]
            labels = sm.slot_labels(f)
            bar_labels = list(labels)
            a = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
            mu, u, s = tn.fresh(tn.LORENTZ), tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
            bar_labels[0] = s
            gen = sm.coupling(r)
            if fac.abelian:
                gen = gen * num(info.charges[r])
            else:
                k = next(i for i, sl in enumerate(info.slots) if sl.factor == r)
                j = tn.fresh(info.slots[k].space)
                bar_labels[k] = j
                gen = gen * tn.T(a[0], j, labels[k])
            Pc, Pcb = projector(sm, f), projector(sm, f, physical=False)
            ref = ref + Density.monomial(F, gen * tn.gammahat(mu, s, u) * Pcb(u, labels[0]), Occurrence(f"gh_{r}", a),
                                         Occurrence(f + "bar", tuple(bar_labels)), Occurrence(f, labels, (mu,)))
            ref = ref + Density.monomial(F, tn.gammahat(mu, s, u) * Pc(u, labels[0]) * gen, Occurrence(f"gh_{r}", a),
                                         Occurrence(f + "bar", tuple(bar_labels), (mu,)), Occurrence(f, labels))
    delta1 = brst.breaking
    assert zero(sm, delta1 - ref)
    assert nonzero(sm, delta1 + ref)
    simp = sm.simplifier()
    assert zero(sm, delta1.physical(simp))
    reduced = LocalFunctional.integrate(ref)
    assert quantum_numbers(sm, reduced) == {(1, Fraction(5))}
    assert sm.rigid_defects(delta1, _field()) == []


@pytest.mark.parametrize("f, r", [("Q", "G"), ("L", "W"), ("E", "B"), ("U", "B")])
def test_breaking_vertex(sm, brst, f, r):
    """report eq:V-delta: V[gh^a(q) Psibar(p') Psi(p)] = -i cal(T)^a (slh(p) P_cbar - slh(p') P_c), p incoming
    (leg 3), p' outgoing (leg 2 has incoming momentum -p')."""
    legs = standard_legs(sm.fields, [f"gh_{r}", f + "bar", f])
    (*a1,), (s2, *g2), (s3, *g3) = (leg.labels for leg in legs)
    value = brst.breaking.vertex_function(legs, sm.simplifier())
    u, nu = tn.fresh(tn.SPINOR), tn.fresh(tn.LORENTZ)
    T = generator(sm, f, r, a1[0] if a1 else None, g2, g3)
    ref = -I * T * tn.gammahat(nu, s2, u) * (tn.p(3, nu) * projector(sm, f, False)(u, s3)
                                             + tn.p(2, nu) * projector(sm, f)(u, s3))
    assert same_expr(sm, value, ref, last_momentum=3)


# ---------------------------------------------------------------------------------------------- identities
def test_linear_identities_and_rho_tilde(sm, brst, simp):
    """report eq:linear-ids and eq:rhotilde: dS_0/db = d.A + xi b; dS_0/d omegabar + d_mu dS_0/d rho^mu = 0; the
    variation of rhotilde = rho + d omegabar has no b (03b: 'the b-terms cancel')."""
    for r in sm.factors:
        assert brst.gauge_condition(brst.action, r).simplify(simp).is_zero(), r
        assert brst.ghost_equation(brst.action, r).simplify(simp).is_zero(), r
        assert brst.ghost_equation(brst.breaking, r).is_zero(), r
        assert not brst.ghost_equation(brst.gauge_fixing, r).simplify(simp).is_zero(), r
        mu, *adj = brst.b_d.images[f"src_{r}"].labels
        b_rho = brst.b_d.images[f"src_{r}"].density
        b_tilde = b_rho + _renamed(brst.b_d.images[f"agh_{r}"], tuple(adj)).total_derivative(mu)
        assert any(f"b_{r}" in c for c in b_rho.field_content()), r
        assert not any(f"b_{r}" in c for c in b_tilde.simplify(simp).field_content()), r


def _renamed(img: FieldImage, labels) -> Density:
    return img.density.rename(list(zip(img.labels, labels)))


def test_rigid_invariance_of_everything(sm, brst):
    """Rigid invariance (report sec:brst): the gauge fermion, L_gf, L_ext, S_0, Delta_1 and the pairings of every
    image of s and b_d with the conjugate field."""
    field = _field()
    for part in (brst.gauge_fermion, brst.gauge_fixing, brst.external, brst.action, brst.breaking):
        assert sm.rigid_defects(part, field) == []
    for name, img in brst.b_d.images.items():
        info = sm.fields[name]
        partner = info.source_of if info.is_source else (info.source or f"b_{info.factor}")
        if info.kind == "antighost":
            continue  # b_d omegabar = b: pairing with b would need a source of omegabar
        assert sm.rigid_defects(paired(sm, partner, img.labels, img.density), field) == [], name


# ---------------------------------------------------------------------------------------------- dictionary, b-table
def test_dictionary_examples(sm, brst):
    """report 03b sec:dictionary examples 1-3 (ghost leg first, all momenta incoming). Examples 1 and 2 hold for b_d
    and b alike; example 3 glues the source vertex to the kinetic vertex, unsplit pslash' with b_d and barred with
    the four-dimensional b."""
    F, b, bd, simp = sm.fields, brst.b, brst.b_d, sm.simplifier()
    gs = sm.coupling("G")
    # 1. Y = -1/2 int (d.A)^2: (i) V(gh(q), A(-q)) = i q^2 q^nu (q^2 four-dimensional, as b is);
    #    (ii) V(gh^c, A^b, A^e) = fg^{bec}[p2 p2 - p1 p1]
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    Y = LocalFunctional.integrate(Density.monomial(F, Expression.num(-1) / 2, Occurrence("B", (mu,), (mu,)),
                                                   Occurrence("B", (nu,), (nu,))))
    legs = standard_legs(F, ["gh_B", "B"])
    assert same_expr(sm, b(Y).vertex_function(legs, simp), I * tn.spbar(1, 1) * tn.p(1, legs[1].labels[0]))
    a = tn.fresh(sm.factors["G"].adjoint_space())
    Y = LocalFunctional.integrate(Density.monomial(F, Expression.num(-1) / 2, Occurrence("G", (mu, a), (mu,)),
                                                   Occurrence("G", (nu, a), (nu,))))
    legs = standard_legs(F, ["gh_G", "G", "G"])
    (c1,), (m2, a2), (m3, a3) = (leg.labels for leg in legs)
    ref = gs * tn.f(a2, a3, c1) * (tn.p(3, m2) * tn.p(3, m3) - tn.p(2, m2) * tn.p(2, m3))
    assert same_expr(sm, bd(Y).vertex_function(legs, simp), ref, last_momentum=3)
    # 2. Y = int psibar gamma^mu Am_mu psi (Q, gluons): V(gh^a(q), Qbar, Q) = -i qslash cal(T)^a P_c
    s, c, w = sm.slot_labels("Q")
    t, u, j = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR), tn.fresh(F["Q"].slots[1].space)
    vertex = LocalFunctional.integrate(Density.monomial(
        F, gs * tn.T(a, c, j) * tn.gammabar(mu, s, u) * tn.PL(u, t), Occurrence("Qbar", (s, c, w)),
        Occurrence("Q", (t, j, w)), Occurrence("G", (mu, a))))
    legs = standard_legs(F, ["gh_G", "Qbar", "Q"])
    (a1,), (s2, c2, w2), (s3, c3, w3) = (leg.labels for leg in legs)
    ref = -I * gs * tn.T(a1, c2, c3) * tn.delta(w2, w3) * tn.gammabar(nu, s2, u) * tn.p(1, nu) * tn.PL(u, s3)
    assert same_expr(sm, bd(vertex).vertex_function(legs, simp), ref)
    assert same_expr(sm, b(vertex).vertex_function(legs, simp), ref)
    # 3. Y = int i Rbar wm psi: V(gh^a, Qbar(p'), Q) = i slb(p') cal(T)^a P_c (p' = -p_2 outgoing),
    #    and with a gluon leg i gamma^mu cal(T)^e cal(T)^a P_c
    Y = LocalFunctional.integrate(Density.monomial(
        F, I * gs * tn.T(a, c, j) * tn.PL(s, t), Occurrence("src_Q", (s, c, w)), Occurrence("gh_G", (a,)),
        Occurrence("Q", (t, j, w))))
    ref = I * gs * tn.T(a1, c2, c3) * tn.delta(w2, w3) * tn.gammabar(nu, s2, u) * (-tn.p(2, nu)) * tn.PL(u, s3)
    assert same_expr(sm, b(Y).vertex_function(legs, simp), ref)
    ref = I * gs * tn.T(a1, c2, c3) * tn.delta(w2, w3) * tn.gamma(nu, s2, u) * (-tn.p(2, nu)) * tn.PL(u, s3)
    assert same_expr(sm, bd(Y).vertex_function(legs, simp), ref)
    legs = standard_legs(F, ["gh_G", "Qbar", "Q", "G"])
    (a1,), (s2, c2, w2), (s3, c3, w3), (m4, a4) = (leg.labels for leg in legs)
    k = tn.fresh(F["Q"].slots[1].space)
    ref = I * gs**2 * tn.gammabar(m4, s2, u) * tn.PL(u, s3) * tn.T(a4, c2, k) * tn.T(a1, k, c3) * tn.delta(w2, w3)
    assert sm.vanishes_modulo_group(simp(b(Y).vertex_function(legs, simp) - ref), [a1, c2, c3, w2, w3, a4], _field())


def test_b_table_entries(sm, brst):
    """report eq:btable-gauge and eq:btable-fermion (four-dimensional b): b Q_M = -N_1, b Q_K^f = sum_r N_F^{fr},
    b Q_V^{fr} = -N_F^{fr} with Q_M = 1/2 A.A, N_1 = int omega d.A, Q_K = int psibar i dslash psi,
    Q_V = int psibar gamma.Am psi, N_F = int psibar gamma^mu (d_mu wm) psi."""
    F, b = sm.fields, brst.b
    for r, fac in sm.factors.items():
        mu = tn.fresh(tn.LORENTZ)
        adj = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        QM = Density.monomial(F, Expression.num(1) / 2, Occurrence(r, (mu,) + adj), Occurrence(r, (mu,) + adj))
        N1 = Density.monomial(F, 1, Occurrence(f"gh_{r}", adj), Occurrence(r, (mu,) + adj, (mu,)))
        assert zero4(sm, b(LocalFunctional.integrate(QM)) + N1), r
    s, c, w = sm.slot_labels("Q")
    t, u, mu = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR), tn.fresh(tn.LORENTZ)
    QK = LocalFunctional.integrate(Density.monomial(F, I * tn.gammabar(mu, s, u) * tn.PL(u, t),
                                                    Occurrence("Qbar", (s, c, w)), Occurrence("Q", (t, c, w), (mu,))))
    NF_sum = Density(F)
    for r in ("G", "W", "B"):
        fac = sm.factors[r]
        a = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        labels_bar, labels = [s, c, w], [t, c, w]
        gen = sm.coupling(r)
        if fac.abelian:
            gen = gen * num(F["Q"].charges[r])
        else:
            k = next(i for i, sl in enumerate(F["Q"].slots) if sl.factor == r)
            j = tn.fresh(F["Q"].slots[k].space)
            labels[k] = j
            gen = gen * tn.T(a[0], labels_bar[k], j)
        struct = gen * tn.gammabar(mu, s, u) * tn.PL(u, t)
        QV = Density.monomial(F, struct, Occurrence("Qbar", tuple(labels_bar)), Occurrence("Q", tuple(labels)),
                              Occurrence(r, (mu,) + a))
        NF = Density.monomial(F, struct, Occurrence("Qbar", tuple(labels_bar)), Occurrence(f"gh_{r}", a, (mu,)),
                              Occurrence("Q", tuple(labels)))
        assert zero4(sm, b(LocalFunctional.integrate(QV)) + NF), r
        NF_sum = NF_sum + NF
    assert zero4(sm, b(QK) - NF_sum)
    assert nonzero4(sm, b(QK) + NF_sum)


# ---------------------------------------------------------------------------------------------- another model
TOY = {"name": "toy", "gauge": {"X": {"group": "SU(2)", "coupling": "gx"}, "Y": {"group": "U(1)", "coupling": "gy"}},
       "fermions": {"P": {"chirality": "L", "X": "F", "Y": "1/3"}, "N": {"chirality": "R", "Y": "1/3"}},
       "scalars": {"phi": {"complex": False, "X": "A"}},
       "parameters": {"mphi": {}},
       "potential": [{"coupling": "mphi", "bilinears": [["phi", "phi"]]}],
       "gauge_fixing": {"kind": "linear_covariant", "xi": "1"}}


def test_a_model_with_a_real_adjoint_scalar():
    """Nothing is specific to the Standard Model: SU(2) x U(1) with a chiral doublet, a singlet and a real adjoint
    scalar (s phi^a = fg^{abc} omega^b phi^c, a real field whose source enters like rho). s^2 = 0, ST(S_0) = s S_0,
    P Delta_1 = 0, S_0 hermitian, b^2 = 0 in four dimensions."""
    from hepkit_renorm import model_from_dict

    m = model_from_dict(TOY)
    toy, simp = BRST(m), m.simplifier()
    assert {p.field: (p.sign, p.side) for p in toy.pairs if p.field in ("phi", "P", "Pbar")} == {
        "phi": (1, "left"), "P": (1, "left"), "Pbar": (-1, "left")}
    S0 = toy.action
    assert zero(m, toy.slavnov_taylor(S0) - toy.breaking)
    assert zero(m, S0 - S0.conjugate())
    assert zero(m, toy.breaking.physical(simp)) and not zero(m, toy.breaking)
    for name, img in toy.s.images.items():
        if m.fields[name].source:
            assert zero(m, paired(m, m.fields[name].source, img.labels, toy.s(img.density))), name
    for name, img in toy.b.images.items():
        info = m.fields[name]
        partner = info.source_of if info.is_source else info.source
        if partner is not None:
            assert zero4(m, paired(m, partner, img.labels, toy.b(img.density))), name
